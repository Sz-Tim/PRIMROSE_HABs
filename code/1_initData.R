# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial dataset compilation


# setup -------------------------------------------------------------------

library(raster)
library(gdistance)
library(tidyverse)
library(glue)
library(lubridate)
library(ncdf4)
library(sf)
library(jsonlite)
library(WeStCOMS)
source("code/00_fn.R")


nDays_avg <- 14
UK_bbox <- list(xmin=-11, xmax=3, ymin=49, ymax=61.5)
sp_i <- read_csv("data/sp_i.csv")
hab.tl <- read_csv("data/hab_tl_thresholds.csv") %>%
  filter(min_ge != -99) %>%
  rename(sp=hab_parameter) %>%
  group_by(sp) %>%
  mutate(alert=case_when(is.na(alert)~0,
                         alert=="warn"~1,
                         alert=="alert"~2),
         A=paste0("A", as.numeric(alert>0)),
         tl=factor(tl, labels=paste0("TL", 0:3))) %>%
  ungroup %>%
  select(sp, min_ge, A, alert, tl)
hab.tl <- hab.tl %>% 
  bind_rows(hab.tl %>% filter(sp=="pseudo_nitzschia_sp") %>% 
              mutate(sp="pseudo_nitzschia_delicatissima_group")) %>% 
  bind_rows(hab.tl %>% filter(sp=="pseudo_nitzschia_sp") %>% 
              mutate(sp="pseudo_nitzschia_seriata_group")) %>%
  mutate(sp=sp_i$abbr[match(sp, sp_i$full)])
# sp_i <- sp_i %>% 
#   left_join(hab.tl %>% filter(A=="A1") %>% 
#               group_by(sp) %>% slice_head(n=1) %>%
#               rename(N_thresh=min_ge, tl_thresh=tl) %>%
#               select(sp, N_thresh, tl_thresh) %>%
#               mutate(N_thresh=pmax(1, N_thresh)) %>%,
#             by=c("abbr"="sp"))
# write_csv(sp_i, "data/sp_i.csv")



# sampling locations and dates --------------------------------------------

fsa.df <- fromJSON("data/0_init/copy_fsa.txt") %>% 
  as_tibble %>% 
  filter(!is.na(date_collected)) %>%
  mutate(datetime_collected=as_datetime(date_collected),
         date=date(datetime_collected)) %>%
  mutate(across(any_of(sp_i$full), ~na_if(.x, -99))) %>%
  filter(datetime_collected >= "2013-07-20") %>%
  group_by(sin) %>%
  mutate(N=n()) %>%
  filter(N > 2) %>%
  ungroup %>%
  mutate(siteid=as.numeric(factor(sin))) %>%
  rename(obsid=oid) %>%
  group_by(sin) %>% 
  mutate(lon=median(easting), lat=median(northing)) %>%
  ungroup %>%
  select(-geom, -easting, -northing, -tide, -datetime_collected, -samplemethod, 
         -depth, -site, -area, -date_collected, -farm_species, -N) %>%
  rename(all_of(setNames(sp_i$full, sp_i$abbr)))

site.df <- fsa.df %>%  
  select(siteid, sin, lon, lat) %>%
  group_by(siteid) %>% slice_head(n=1) %>% ungroup
# saveRDS(site.df, "data/site_df.rds")

fsa.df <- fsa.df %>% select(-lon, -lat, -sin)
saveRDS(fsa.df, "data/0_init/fsa_df.rds")




# CMEMS -------------------------------------------------------------------

# Atlantic-European North West Shelf Biogeochemistry Analysis and Forecast
# Reanalysis: 1993-present, but updated every 6 months
#  - https://doi.org/10.48670/moi-00058
# Analysis & Forecast: 2019-present, updated daily with 6-day forecast
#  - https://doi.org/10.48670/moi-00056

cmems_i <- expand_grid(
  var=c("chl", # Mass concentration of chlorophyll a
        "kd", # Volume beam attenuation coefficient of radiative flux
        "no3", # Mole concentration of nitrate
        "o2", # Mole concentration of dissolved molecular oxygen
        "ph", # Sea water ph reported on total scale
        "phyc", # Mole concentration of phytoplankton expressed as carbon
        "po4", # Mole concentration of phosphate
        "pp", # Net primary production of biomass expressed as carbon per unit volume
        "spco2" # Surface partial pressure of carbon dioxide
  ),
  source=c("Reanalysis", "Analysis&Forecast")) %>%
  mutate(server=if_else(source=="Reanalysis", "@my.cmems-du.eu", "@nrt.cmems-du.eu"), 
         doi=glue("https://doi.org/10.48670/moi-0005{if_else(source=='Reanalysis', 8, 6)}"),
         ID=glue("cmems_mod_nws_bgc-{var}_", 
                 "{if_else(source=='Reanalysis', 'my', 'anfc')}_7km-",
                 "{if_else(var=='spco2', '2D', '3D')}_P1D-m"))
write_csv(cmems_i, "data/cmems_i.csv")

cmems_cred <- readRDS("data/cmems_cred.rds")
get_CMEMS(userid=cmems_cred$userid, pw=cmems_cred$pw, 
          i.df=cmems_i, bbox=UK_bbox, 
          nDays_buffer=nDays_avg, dateRng=range(fsa.df$date), 
          out.dir="data/0_init/cmems/")

rean.f <- dirf("data/0_init/cmems", "cmems.*Reanalysis.rds")
anfo.f <- dirf("data/0_init/cmems", "cmems.*Forecast.rds")
cmems.df <- bind_rows(
  readRDS(rean.f[1]) %>%
    bind_cols(map_dfc(rean.f[-1], ~readRDS(.x) %>% select(5))) %>%
    filter(date < "2019-05-01"),
  readRDS(anfo.f[1]) %>%
    bind_cols(map_dfc(anfo.f[-1], ~readRDS(.x) %>% select(5)))
  ) %>%
  na.omit() %>%
  mutate(chl=log1p(chl),
         kd=log(kd),
         no3=log1p(no3),
         o2=log(o2),
         phyc=log1p(phyc),
         po4=log1p(po4)) %>%
  arrange(date, lon, lat) %>%
  group_by(date) %>%
  mutate(cmems_id=row_number()) %>%
  ungroup
saveRDS(cmems.df, glue("data/0_init/cmems_end_{max(cmems.df$date)}.rds"))




# WeStCOMS-WRF ------------------------------------------------------------

# There are three different resolutions / domains, from coarse but expansive
# (d01) to finer but restricted (d03). The getWRF() function has an argument to
# nest all three, or select just one. WeStCOMS-WRF
#  - Wind speed
#  - Wind direction
#  - Shortwave radiation
#  - Precipitation
#  - Sea surface temperature

wrf.dir <- ifelse(.Platform$OS.type=="unix",
                  "/media/archiver/common/sa01da-work/WRF/Archive/",
                  "D:/hydroOut/WRF/Archive/")
get_WRF(wrf.dir=wrf.dir, nDays_buffer=nDays_avg, 
        dateRng=c(ymd("2016-01-07"), max(fsa.df$date)), 
        out.dir="data/0_init/")

# read and subset WRF domains to nest higher res within lower res
wrf.out <- "data/0_init/wrf/"
wrf.df <- subset_WRF("d01", wrf.out, v2_start="2019-04-01") %>%
  bind_rows(subset_WRF("d02", wrf.out, v2_start="2019-04-01")) %>%
  bind_rows(subset_WRF("d03", wrf.out, v2_start="2019-04-01")) %>%
  arrange(date, res, i) %>%
  group_by(date) %>%
  mutate(wrf_id=row_number()) %>%
  ungroup %>%
  mutate(Shortwave=log1p(Shortwave),
         UV_mn=log1p(UV_mn),
         Precip=log1p(pmax(Precip, 0)*3600*24*1000), # m/s to mm/day
         sst=if_else(sst > -100, sst, NA_real_))
saveRDS(wrf.df, glue("data/0_init/wrf_end_{max(wrf.df$date)}.rds"))




# Pairwise distances ------------------------------------------------------

# Shortest paths within the ocean

path.ls <- get_shortestPaths(ocean.path="data/northAtlantic_footprint.tif", 
                             site.df=site.df, 
                             transMx.path="data/mesh_tmx.rds", recalc_transMx=F)
write_csv(path.ls$dist.df, "data/site_pairwise_distances.csv")
site.df <- path.ls$site.df # slightly modified lat/lon to fit within ocean mesh
path.ls$dist.df %>% 
  bind_rows(tibble(origin=unique(.$origin), 
                   destination=unique(.$origin), 
                   distance=0)) %>%
  filter(distance < 100e3) %>% 
  select(-distance) %>% 
  group_by(origin) %>% 
  nest(data=destination) %>%
  mutate(dest_c=c(data[[1]])) %>% 
  select(-data) %>%
  ungroup %>%
  saveRDS("data/site_neighbors_100km.rds")




# Fetch and bearing -------------------------------------------------------

# Wave fetch and bearing with the most open water
# https://doi.org/10.6084/m9.figshare.12029682.v1

site.df <- site.df %>%
  get_fetch(., "data/log10_eu200m1a.tif") %>%
  get_openBearing(., "data/northAtlantic_footprint.gpkg", buffer=200e3)
# saveRDS(site.df, "data/site_df.rds")





# site buffers ------------------------------------------------------------

# Buffers for averaging environmental conditions
# CMEMS is coarse and does not resolve lochs, so some sites are not covered
# The buffer size was determined by trial and error to reduce averaging while
# ensuring all sites are represented.

site.df <- readRDS("data/site_df.rds")
site.sf <- site.df %>% 
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  select(-sin) %>%
  st_buffer(dist=100e3) %>%
  split_to_NSEW()
st_write(site.sf, "data/site_sf.gpkg")





# HAB densities -----------------------------------------------------------

# Observed densities and bloom states for focal taxa
# * = [t]
# *1 = [t-1]
# *2 = [t-2]
# N: reported density
# lnN: ln(N + 1)
# tl: HABReports traffic light for density N
# cat: HABReports label category for density N
# alert: HABReports action (0_none, 1_warn, 2_alert)
# lnNAvg: regional average of lnN within previous week
# prAlertAvg: regional average of alerting sites within previous week

site.df <- readRDS("data/site_df.rds")
hab.df <- calc_hab_features(readRDS("data/0_init/fsa_df.rds"),
                            sp_i, hab.tl,
                            readRDS("data/site_neighbors_100km.rds"))
saveRDS(hab.df, "data/0_init/hab_densities.rds")





# Extract sites -----------------------------------------------------------

# . CMEMS  site:date ------------------------------------------------------

site.df <- readRDS("data/site_df.rds")
cmems_vars <- c("chl", "kd", "no3", "o2", "ph", "phyc", "po4", "pp")
cmems.df <- readRDS(dirf("data/0_init", "cmems_end.*rds"))
cmems.sf <- cmems.df %>% select(date, lon, lat, cmems_id) %>% 
  filter(date==min(date)) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326)
site.df <- site.df %>%
  st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) %>%
  st_transform(4326) %>%
  mutate(cmems_id=st_nearest_feature(., cmems.sf)) %>%
  st_drop_geometry()
saveRDS(site.df, "data/site_df.rds")
cmems.site <- cmems.df %>% 
  filter(cmems_id %in% site.df$cmems_id) %>%
  group_by(cmems_id) %>%
  mutate(across(all_of(cmems_vars), 
                ~zoo::rollmean(.x, k=7, na.pad=T),
                .names="{.col}Wk")) %>%
  mutate(across(any_of(paste0(cmems_vars, "Wk")),
                ~.x - lag(.x),
                .names="{.col}Delta")) %>%
  ungroup %>%
  mutate(yday=yday(date)) %>%
  group_by(cmems_id) %>%
  mutate(across(all_of(cmems_vars),
                ~detrend_loess(yday, .x, span=0.3), 
                .names="{.col}Dt")) %>%
  ungroup
cmems.site <- cmems.site %>% 
  select(cmems_id, date, 
         all_of(paste0(cmems_vars, "Wk")),
         all_of(paste0(cmems_vars, "WkDelta")),
         all_of(paste0(cmems_vars, "Dt")))
saveRDS(cmems.site, "data/0_init/cmems_sitePt.rds")


# . CMEMS  buffer:date ----------------------------------------------------

site.buffer <- st_read("data/site_sf.gpkg") %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  mutate(cmems_id=st_intersects(., cmems.sf)) %>%
  st_drop_geometry()
cmems.buffer <- expand_grid(siteid=site.df$siteid,
                            quadrant=unique(site.buffer$quadrant),
                            date=unique(cmems.df$date)) %>%
  bind_cols(as_tibble(setNames(map(cmems_vars, ~NA_real_), cmems_vars)))

cmems_id.ls <- map(site.buffer$cmems_id, ~.x)
cmems_dates.ls <- map(unique(cmems.df$date), ~which(cmems.df$date==.x))
ij <- 1
for(i in 1:nrow(site.buffer)) {
  for(j in 1:length(cmems_dates.ls)) {
    if(length(cmems_id.ls[[i]]) > 0) {
      for(k in cmems_vars) {
        cmems.buffer[ij,k] <- mean(cmems.df[cmems_dates.ls[[j]],][cmems_id.ls[[i]],][[k]])
      } 
    }
    ij <- ij+1
    if(ij %% 1000 == 0) {cat(ij, "of", nrow(cmems.buffer), "\n")}
  }
}
cmems.buffer <- cmems.buffer %>% 
  group_by(siteid, quadrant) %>%
  mutate(across(any_of(cmems_vars), 
                ~zoo::rollmean(.x, k=7, na.pad=T),
                .names="{.col}AvgWk")) %>%
  mutate(across(any_of(paste0(cmems_vars, "AvgWk")),
                ~.x - lag(.x),
                .names="{.col}Delta")) %>%
  ungroup %>%
  mutate(yday=yday(date)) %>%
  group_by(siteid, quadrant) %>%
  mutate(across(any_of(cmems_vars),
                ~detrend_loess(yday, .x, span=0.3), 
                .names="{.col}AvgDt")) %>%
  ungroup %>% 
  select(siteid, quadrant, date, 
         all_of(paste0(cmems_vars, "AvgWk")),
         all_of(paste0(cmems_vars, "AvgWkDelta")),
         all_of(paste0(cmems_vars, "AvgDt"))) %>%
  group_by(siteid, date) %>%
  mutate(across(where(is.numeric), zoo::na.aggregate)) %>%
  ungroup
saveRDS(cmems.buffer, "data/0_init/cmems_siteBufferNSEW.rds")



# . WRF  site:date --------------------------------------------------------

wrf_i <- list(vars=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea_vars=c("U", "V", "UV", "Shortwave", "Precip"),
              land_vars=c("sst"))
site.df <- readRDS("data/site_df.rds") %>% select(-starts_with("wrf_id"))
wrf_versions <- map(seq_along(dir("data/0_init/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/0_init/wrf", glue("domain_d0._{.x}")), readRDS) %>%
                      arrange(res, i) %>%
                      mutate(wrf_id=row_number()) %>%
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
site.df <- map(wrf_versions,
               ~site.df %>%
                 st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) %>%
                 st_transform(4326) %>%
                 mutate(wrf_id=st_nearest_feature(., .x)) %>%
                 st_drop_geometry) %>%
  reduce(full_join, 
         by=names(site.df), suffix=paste0(".", seq_along(wrf_versions)))
# saveRDS(site.df, "data/site_df.rds")
site.versions <- grep("wrf_id", names(site.df), value=T)
wrf.df <- readRDS(dirf("data/0_init/", "wrf_end_.*rds")) %>%
  mutate(sst=if_else(sst > -100, sst, NA_real_))
wrf.site <- wrf.df %>%
  group_by(version) %>%
  group_split() %>%
  map2_dfr(., site.versions, 
           ~.x %>% filter(wrf_id %in% site.df[[.y]])) %>%
  arrange(wrf_id, date) %>%
  group_by(wrf_id) %>%
  mutate(across(any_of(wrf_i$var), 
                ~zoo::rollmean(.x, k=7, na.pad=T),
                .names="{.col}Wk")) %>%
  mutate(across(any_of(paste0(wrf_i$var, "Wk")),
                ~.x - lag(.x),
                .names="{.col}Delta")) %>%
  ungroup %>%
  mutate(yday=yday(date)) %>%
  group_by(wrf_id) %>%
  mutate(across(any_of(wrf_i$var),
                ~detrend_loess(yday, .x, span=0.3), 
                .names="{.col}Dt")) %>%
  ungroup
wrf.site <- wrf.site %>% 
  select(wrf_id, version, date, 
         all_of(paste0(wrf_i$vars, "Wk")),
         all_of(paste0(wrf_i$vars, "WkDelta")),
         all_of(paste0(wrf_i$vars, "Dt")))
saveRDS(wrf.site, "data/0_init/wrf_sitePt.rds")




# . WRF  buffer:date ------------------------------------------------------

site.buffer <- map(wrf_versions,
                   ~st_read("data/site_sf.gpkg") %>%
                     select(siteid, quadrant, geom) %>%
                     st_transform(4326) %>%
                     st_make_valid() %>%
                     mutate(wrf_id=st_intersects(., .x)) %>%
                     st_drop_geometry()) %>%
  reduce(full_join, 
         by=c("siteid", "quadrant"), 
         suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer <- expand_grid(siteid=site.df$siteid,
                          quadrant=unique(site.buffer$quadrant),
                          date=unique(wrf.df$date)) %>%
  mutate(version=1 + (date >= "2019-04-01")) %>%
  bind_cols(as_tibble(setNames(map(wrf_i$vars, ~NA_real_), wrf_i$vars)))

wrf_id.ls <- list(v1=map(site.buffer$wrf_id.1, ~.x),
                  v2=map(site.buffer$wrf_id.2, ~.x))
wrf_dates.ls <- map(unique(wrf.df$date), ~which(wrf.df$date==.x))
ij <- 1
for(i in 1:nrow(site.buffer)) {
  for(j in 1:length(wrf_dates.ls)) {
    if(length(wrf_id.ls[[wrf.buffer$version[ij]]][[i]]) > 0) {
      for(k in wrf_i$sea_vars) {
        wrf.buffer[ij,k] <- mean(wrf.df[wrf_dates.ls[[j]],][wrf_id.ls[[wrf.buffer$version[ij]]][[i]],][[k]])
      }
      sst_ij <- wrf.df[wrf_dates.ls[[j]],][wrf_id.ls[[wrf.buffer$version[ij]]][[i]],][["sst"]]
      elev_ij <- wrf.df[wrf_dates.ls[[j]],][wrf_id.ls[[wrf.buffer$version[ij]]][[i]],][["elev"]]
      wrf.buffer[ij,"sst"] <- mean(sst_ij[elev_ij==0], na.rm=T) 
    }
    ij <- ij+1
    if(ij %% 1000 == 0) {cat(ij, "of", nrow(wrf.buffer), "\n")}
  }
}
wrf.buffer <- wrf.buffer %>% 
  group_by(siteid, quadrant) %>%
  mutate(across(any_of(wrf_i$vars), 
                ~zoo::rollmean(.x, k=7, na.pad=T),
                .names="{.col}AvgWk")) %>%
  mutate(across(any_of(paste0(wrf_i$vars, "AvgWk")),
                ~.x - lag(.x),
                .names="{.col}Delta")) %>%
  ungroup %>%
  mutate(yday=yday(date)) %>%
  group_by(siteid, quadrant) %>%
  mutate(across(any_of(wrf_i$vars),
                ~detrend_loess(yday, .x, span=0.3), 
                .names="{.col}AvgDt")) %>%
  ungroup %>% 
  select(siteid, quadrant, date, 
         all_of(paste0(wrf_i$vars, "AvgWk")),
         all_of(paste0(wrf_i$vars, "AvgWkDelta")),
         all_of(paste0(wrf_i$vars, "AvgDt"))) %>%
  group_by(siteid, date) %>%
  mutate(across(where(is.numeric), zoo::na.aggregate)) %>%
  ungroup
saveRDS(wrf.buffer, "data/0_init/wrf_siteBufferNSEW.rds")




# Compile -----------------------------------------------------------------

# Load all datasets and merge
site.df <- readRDS("data/site_df.rds") %>% select(-sin)
hab.df <- readRDS("data/0_init/hab_densities.rds")
cmems.site <- readRDS("data/0_init/cmems_sitePt.rds")
cmems.buffer <- readRDS("data/0_init/cmems_siteBufferNSEW.rds") %>% 
  pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")
wrf.site <- readRDS("data/0_init/wrf_sitePt.rds") %>% select(-version, -contains("sst"))
wrf.buffer <- readRDS("data/0_init/wrf_siteBufferNSEW.rds") %>% 
  pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")

# load, extract, and compile
# - Algal densities
# - CMEMS
# - WRF
# - fetch & bearing
# - calculate cos(dir)'s
obs.df <- site.df %>%
  right_join(hab.df, by="siteid", multiple="all") %>%
  left_join(cmems.site, by=c("cmems_id", "date")) %>%
  left_join(cmems.buffer, by=c("siteid", "date")) %>%
  mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) %>%
  select(-wrf_id.1, wrf_id.2) %>%
  left_join(wrf.site, by=c("wrf_id", "date")) %>%
  left_join(wrf.buffer, by=c("siteid", "date")) %>%
  na.omit()
saveRDS(obs.df, "data/0_init/data_full_allSpp.rds")

cmems_vars <- grep("cmems_id|date|siteid", 
                   unique(c(names(cmems.site), names(cmems.buffer))),
                   value=T, invert=T) %>% sort
saveRDS(cmems_vars, "data/cmems_vars.rds")

wrf_vars <- grep("wrf_id|date|siteid", 
                   unique(c(names(wrf.site), names(wrf.buffer))),
                   value=T, invert=T) %>% sort
saveRDS(wrf_vars, "data/wrf_vars.rds")


# Reduce highly correlated predictors
corrplot::corrplot(cor(obs.df[,wrf_vars], use="pairwise"), diag=F, method="number", 
                   order="alphabet", number.cex=0.7, tl.cex=0.6)
corrplot::corrplot(abs(cor(obs.df[,wrf_vars], use="pairwise"))>0.9, is.corr=F, diag=F,
                   order="alphabet", number.cex=0.7, tl.cex=0.8)
corrplot::corrplot(cor(obs.df[,cmems_vars], use="pairwise"), diag=F, method="number", 
                   order="alphabet", number.cex=0.7, tl.cex=0.6)
corrplot::corrplot(abs(cor(obs.df[,cmems_vars], use="pairwise"))>0.9, is.corr=F, diag=F,
                   order="alphabet", number.cex=0.7, tl.cex=0.8)




