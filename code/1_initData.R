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
tox_i <- read_csv("data/tox_i.csv")
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
# TODO: what do they use for the alerts?
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

fsa.df <- fromJSON("data/0_init/copy_fsa.txt") %>% as_tibble %>% 
  filter(!is.na(date_collected)) %>%
  mutate(datetime_collected=as_datetime(date_collected),
         date=date(datetime_collected)) %>%
  mutate(across(any_of(sp_i$full), ~na_if(.x, -99))) %>%
  filter(datetime_collected >= "2013-07-20") %>%
  group_by(sin) %>% mutate(N=n()) %>% filter(N > 2) %>% ungroup %>%
  mutate(siteid=as.numeric(factor(sin))) %>%
  rename(obsid=oid) %>%
  group_by(sin) %>% mutate(lon=median(easting), lat=median(northing)) %>% ungroup %>%
  select(-geom, -easting, -northing, -tide, -datetime_collected, -samplemethod, 
         -depth, -site, -area, -date_collected, -farm_species, -N) %>%
  rename(all_of(setNames(sp_i$full, sp_i$abbr)))
site_hab.df <- fsa.df %>%  
  select(siteid, sin, lon, lat) %>%
  group_by(siteid) %>% slice_head(n=1) %>% ungroup
# saveRDS(site_hab.df, "data/site_hab_df.rds")
fsa.df <- fsa.df %>% select(-lon, -lat, -sin)
saveRDS(fsa.df, "data/0_init/fsa_df.rds")

# TODO: What are the negative numbers? What are the alert thresholds?
cefas.df <- fromJSON("data/0_init/copy_cefas.txt") %>% as_tibble %>%
  filter(sin != "-99", !is.na(date_collected)) %>%
  mutate(datetime_collected=as_datetime(date_collected),
         date=date(datetime_collected),
         across(one_of("psp", "oa_dtxs_ptxs", "azas", "ytxs", "asp"),
                ~if_else(.x < 0, NA_real_, .x))) %>%
  group_by(sin) %>% mutate(N=n()) %>% filter(N > 2) %>% ungroup %>%
  mutate(siteid=as.numeric(factor(sin))) %>%
  rename(obsid=oid) %>%
  group_by(sin) %>% mutate(lon=median(easting), lat=median(northing)) %>% ungroup %>%
  filter(lat > 500000) %>% # Scotland 
  select(-geom, -easting, -northing, -datetime_collected, 
         -site, -area, -date_collected, -farm_species, -N)
site_tox.df <- cefas.df %>%  
  select(siteid, sin, lon, lat) %>%
  group_by(siteid) %>% slice_head(n=1) %>% ungroup
saveRDS(site_tox.df, "data/site_tox_df.rds")
cefas.df <- cefas.df %>% select(-lon, -lat, -sin)
saveRDS(cefas.df, "data/0_init/cefas_df.rds")



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
                             site.df=site_hab.df, 
                             transMx.path="data/mesh_tmx.rds", recalc_transMx=F)
write_csv(path.ls$dist.df, "data/site_hab_pairwise_distances.csv")
site_hab.df <- path.ls$site.df # slightly modified lat/lon to fit within ocean mesh
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
  saveRDS("data/site_hab_neighbors_100km.rds")

path.ls <- get_shortestPaths(ocean.path="data/northAtlantic_footprint.tif", 
                             site.df=site_tox.df, 
                             transMx.path="data/mesh_tmx.rds", recalc_transMx=F)
write_csv(path.ls$dist.df, "data/site_tox_pairwise_distances.csv")
site_tox.df <- path.ls$site.df # slightly modified lat/lon to fit within ocean mesh
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
  saveRDS("data/site_tox_neighbors_100km.rds")



# Fetch and bearing -------------------------------------------------------

# Wave fetch and bearing with the most open water
# https://doi.org/10.6084/m9.figshare.12029682.v1

site_hab.df <- site_hab.df %>%
  get_fetch(., "data/log10_eu200m1a.tif") %>%
  get_openBearing(., "data/northAtlantic_footprint.gpkg", buffer=200e3)
# saveRDS(site_hab.df, "data/site_hab_df.rds")

site_tox.df <- site_tox.df %>%
  get_fetch(., "data/log10_eu200m1a.tif") %>%
  get_openBearing(., "data/northAtlantic_footprint.gpkg", buffer=200e3)
# saveRDS(site_tox.df, "data/site_tox_df.rds")



# site buffers ------------------------------------------------------------

# Buffers for averaging environmental conditions
# CMEMS is coarse and does not resolve lochs, so some sites are not covered
# The buffer size was determined by trial and error to reduce averaging while
# ensuring all sites are represented.

site_hab.df <- readRDS("data/site_hab_df.rds")
site_hab.sf <- site_hab.df %>% 
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  select(-sin) %>%
  st_buffer(dist=100e3) %>%
  split_to_NSEW()
st_write(site_hab.sf, "data/site_hab_sf.gpkg")

site_tox.df <- readRDS("data/site_tox_df.rds")
site_tox.sf <- site_tox.df %>% 
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  select(-sin) %>%
  st_buffer(dist=100e3) %>%
  split_to_NSEW()
st_write(site_tox.sf, "data/site_tox_sf.gpkg")



# HAB densities -----------------------------------------------------------

# Observed densities and bloom states for focal taxa
# * = [t]
# *1 = [t-1]
# *2 = [t-2]
# N: reported density
# lnN: ln(N + 1)
# tl: HABReports traffic light for density N
# alert: HABReports action (0_none, 1_warn, 2_alert)
# lnNAvg: regional average of lnN within previous week
# prAlertAvg: regional average of alerting sites within previous week

hab.df <- calc_hab_features(readRDS("data/0_init/fsa_df.rds"),
                            sp_i, hab.tl,
                            readRDS("data/site_hab_neighbors_100km.rds"))
saveRDS(hab.df, "data/0_init/hab_obs.rds")



# toxin concentrations ----------------------------------------------------

# Observed concentrations and alert states for focal toxins
# * = [t]
# *1 = [t-1]
# *2 = [t-2]
# N: reported concentration
# lnN: ln(N + 1)
# tl: HABReports traffic light for concentration N
# alert: HABReports action (0_none, 1_warn, 2_alert)
# lnNAvg: regional average of lnN within previous week
# prAlertAvg: regional average of alerting sites within previous week

tox.df <- calc_toxin_features(readRDS("data/0_init/cefas_df.rds"),
                              sp_i, tox.tl,
                              readRDS("data/site_tox_neighbors_100km.rds"))
saveRDS(tox.df, "data/0_init/tox_obs.rds")



# Extract sites -----------------------------------------------------------

# . CMEMS  site:date ------------------------------------------------------
cmems_vars <- c("chl", "kd", "no3", "o2", "ph", "phyc", "po4", "pp")
cmems.df <- readRDS(dirf("data/0_init", "cmems_end.*rds"))
cmems.sf <- cmems.df %>% select(date, lon, lat, cmems_id) %>% 
  filter(date==min(date)) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326)

# HABs
site_hab.df <- readRDS("data/site_hab_df.rds")
site_hab.df <- site_hab.df %>% find_nearest_feature_id(cmems.sf, "cmems_id")
saveRDS(site_hab.df, "data/site_hab_df.rds")
cmems.site_hab <- extract_cmems_pts(site_hab.df, cmems_vars, cmems.df)
saveRDS(cmems.site_hab, "data/0_init/cmems_sitePt_hab.rds")

# toxins
site_tox.df <- readRDS("data/site_tox_df.rds")
site_tox.df <- site_tox.df %>% find_nearest_feature_id(cmems.sf, "cmems_id")
saveRDS(site_tox.df, "data/site_tox_df.rds")
cmems.site_tox <- extract_cmems_pts(site_tox.df, cmems_vars, cmems.df)
saveRDS(cmems.site_hab, "data/0_init/cmems_sitePt_tox.rds")



# . CMEMS  buffer:date ----------------------------------------------------
cmems_vars <- c("chl", "kd", "no3", "o2", "ph", "phyc", "po4", "pp")
cmems.df <- readRDS(dirf("data/0_init", "cmems_end.*rds"))
cmems.sf <- cmems.df %>% select(date, lon, lat, cmems_id) %>% 
  filter(date==min(date)) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326)

# HABs
site.buffer_hab <- st_read("data/site_hab_sf.gpkg") %>%
  find_buffer_intersect_ids(., cmems.sf, "cmems_id")
cmems.buffer_hab <- extract_cmems_buffers(site.buffer_hab, cmems_vars, cmems.df)
saveRDS(cmems.buffer_hab, "data/0_init/cmems_siteBufferNSEW_hab.rds")

# toxins
site.buffer_tox <- st_read("data/site_tox_sf.gpkg") %>%
  find_buffer_intersect_ids(., cmems.sf, "cmems_id")
cmems.buffer_tox <- extract_cmems_buffers(site.buffer_tox, cmems_vars, cmems.df)
saveRDS(cmems.buffer_tox, "data/0_init/cmems_siteBufferNSEW_tox.rds")




# . WRF  site:date --------------------------------------------------------
wrf_i <- list(vars=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea_vars=c("U", "V", "UV", "Shortwave", "Precip"),
              land_vars=c("sst"))
wrf_versions <- map(seq_along(dir("data/0_init/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/0_init/wrf", glue("domain_d0._{.x}")), readRDS) %>%
                      arrange(res, i) %>%
                      mutate(wrf_id=row_number()) %>%
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(dirf("data/0_init/", "wrf_end_.*rds")) %>%
  mutate(sst=if_else(sst > -100, sst, NA_real_))

# HABs
site_hab.df <- readRDS("data/site_hab_df.rds") %>% select(-starts_with("wrf_id"))
site_hab.df <- map(wrf_versions, ~site_hab.df %>% find_nearest_feature(., .x, "wrf_id")) %>%
  reduce(full_join, by=names(site_hab.df), suffix=paste0(".", seq_along(wrf_versions)))
saveRDS(site_hab.df, "data/site_hab_df.rds")
site_hab.versions <- grep("wrf_id", names(site_hab.df), value=T)
wrf.site_hab <- extract_wrf_pts(site_hab.df, site_hab.versions, wrf_i, wrf.df)
saveRDS(wrf.site_hab, "data/0_init/wrf_sitePt_hab.rds")

# toxins
site_tox.df <- readRDS("data/site_tox_df.rds") %>% select(-starts_with("wrf_id"))
site_tox.df <- map(wrf_versions, ~site_tox.df %>% find_nearest_feature(., .x, "wrf_id")) %>%
  reduce(full_join, by=names(site_tox.df), suffix=paste0(".", seq_along(wrf_versions)))
saveRDS(site_tox.df, "data/site_tox_df.rds")
site_tox.versions <- grep("wrf_id", names(site_tox.df), value=T)
wrf.site_tox <- extract_wrf_pts(site_tox.df, site_tox.versions, wrf_i, wrf.df)
saveRDS(wrf.site_tox, "data/0_init/wrf_sitePt_tox.rds")



# . WRF  buffer:date ------------------------------------------------------
wrf_i <- list(vars=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea_vars=c("U", "V", "UV", "Shortwave", "Precip"),
              land_vars=c("sst"))
wrf_versions <- map(seq_along(dir("data/0_init/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/0_init/wrf", glue("domain_d0._{.x}")), readRDS) %>%
                      arrange(res, i) %>%
                      mutate(wrf_id=row_number()) %>%
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(dirf("data/0_init/", "wrf_end_.*rds")) %>%
  mutate(sst=if_else(sst > -100, sst, NA_real_))

# HABs
site.buffer_hab <- map(wrf_versions, 
                       ~st_read("data/site_hab_sf.gpkg") %>% 
                         find_buffer_intersect_ids(., .x, "wrf_id")) %>%
  reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer_hab <- extract_wrf_buffers(site.buffer_hab, wrf_i, wrf.df)
saveRDS(wrf.buffer_hab, "data/0_init/wrf_siteBufferNSEW_hab.rds")

# toxins
site.buffer_tox <- map(wrf_versions, 
                       ~st_read("data/site_tox_sf.gpkg") %>% 
                         find_buffer_intersect_ids(., .x, "wrf_id")) %>%
  reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer_tox <- extract_wrf_buffers(site.buffer_tox, wrf_i, wrf.df)
saveRDS(wrf.buffer_tox, "data/0_init/wrf_siteBufferNSEW_tox.rds")



# Compile -----------------------------------------------------------------

# HABs
hab.ls <- load_datasets("0_init", "hab")
hab.df <- hab.ls$site %>% select(-sin) %>%
  right_join(hab.ls$obs, by="siteid", multiple="all") %>%
  left_join(hab.ls$cmems.pt, by=c("cmems_id", "date")) %>%
  left_join(hab.ls$cmems.buf %>% 
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) %>%
  select(-wrf_id.1, wrf_id.2) %>%
  left_join(hab.ls$wrf.pt, by=c("wrf_id", "date")) %>%
  left_join(hab.ls$wrf.buf %>% 
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  na.omit()
saveRDS(hab.df, "data/0_init/data_hab_all.rds")

# toxins
tox.ls <- load_datasets("0_init", "tox")
tox.df <- tox.ls$site %>% select(-sin) %>%
  right_join(tox.ls$obs, by="siteid", multiple="all") %>%
  left_join(tox.ls$cmems.pt, by=c("cmems_id", "date")) %>%
  left_join(tox.ls$cmems.buf %>% 
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) %>%
  select(-wrf_id.1, wrf_id.2) %>%
  left_join(tox.ls$wrf.pt, by=c("wrf_id", "date")) %>%
  left_join(tox.ls$wrf.buf %>% 
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  na.omit()
saveRDS(tox.df, "data/0_init/data_tox_all.rds")

# variable names
grep("cmems_id|date|siteid", 
     unique(c(names(hab.ls$cmems.pt), names(hab.ls$cmems.buf))),
     value=T, invert=T) %>% 
  sort %>%
  saveRDS("data/cmems_vars.rds")
grep("wrf_id|date|siteid", 
     unique(c(names(hab.ls$wrf.pt), names(hab.ls$wrf.buf))),
     value=T, invert=T) %>% 
  sort %>%
  saveRDS("data/wrf_vars.rds")


