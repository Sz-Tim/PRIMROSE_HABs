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
dateStart <- "2016-01-01"
UK_bbox <- list(xmin=-11, xmax=3, ymin=49, ymax=61.5)
urls <- c(fsa="fsa_counts", 
          cefas="cefas_counts",
          sites="all_sites") %>%
  map(~glue("http://varro:3001/{.x}"))

hab_i <- read_csv("data/i_hab.csv")
tox_i <- read_csv("data/i_tox.csv")
tl_hab <- read_csv("data/tl_thresholds_hab.csv") %>%
  filter(min_ge != -99) %>%
  group_by(abbr) %>%
  mutate(alert=case_when(is.na(alert)~0,
                         alert=="warn"~1,
                         alert=="alert"~2),
         A=paste0("A", as.numeric(alert>0)),
         tl=factor(tl)) %>%
  ungroup %>%
  select(abbr, min_ge, A, alert, tl)
# hab_i <- hab_i %>% select(-ends_with("thresh")) %>%
#   left_join(tl_hab %>% filter(A=="A1") %>%
#               group_by(abbr) %>% slice_head(n=1) %>%
#               rename(N_thresh=min_ge, tl_thresh=tl) %>%
#               select(abbr, N_thresh, tl_thresh))# %>%
# # mutate(N_thresh=pmax(1, N_thresh))) # Why did I do this?
# write_csv(hab_i, "data/i_hab.csv")
# NB: toxin thresholds are adjusted for ASP, AZAs, YTXs because of the extremely 
# limited occurrence of TL2 and TL3. Models thus predict only presence/absence
tl_tox <- read_csv("data/tl_thresholds_tox.csv") %>%
  filter(min_ge != -99) %>%
  group_by(abbr) %>%
  mutate(alert=case_when(is.na(alert)~0,
                         alert=="warn"~1,
                         alert=="alert"~2),
         A=paste0("A", as.numeric(alert>0)),
         tl=factor(tl)) %>%
  ungroup %>%
  select(abbr, min_ge, A, alert, tl) %>%
  mutate(A=if_else(abbr %in% c("asp", "azas", "ytxs") & as.numeric(tl)>1, "A1", A),
         alert=if_else(abbr %in% c("asp", "azas", "ytxs") & as.numeric(tl)>1, 2, alert))
# tox_i <- tox_i %>% select(-ends_with("thresh")) %>%
#   left_join(tl_tox %>% filter(A=="A1") %>%
#               group_by(abbr) %>% slice_head(n=1) %>%
#               rename(N_thresh=min_ge, tl_thresh=tl) %>%
#               select(abbr, N_thresh, tl_thresh))# %>%
#               # mutate(N_thresh=pmax(1, N_thresh))) # Why did I do this?
# write_csv(tox_i, "data/i_tox.csv")



# sampling locations and dates --------------------------------------------

sites <- fromJSON(readLines(url(urls$sites), warn=F)) %>% as_tibble %>%
  filter(east < 2e6) %>%
  mutate(fromdate=date(fromdate), todate=date(todate)) %>%
  rowwise() %>%
  mutate(date=list(seq(fromdate, todate, by=1))) %>%
  ungroup %>%
  select(sin, east, north, date) %>%
  unnest(date)

fsa.df <- fromJSON(readLines(url(urls$fsa), warn=F)) %>% as_tibble %>% 
  filter(!is.na(date_collected)) %>%
  mutate(datetime_collected=as_datetime(date_collected),
         date=date(datetime_collected)) %>%
  filter(date > dateStart) %>%
  mutate(across(any_of(hab_i$full), ~na_if(.x, -99))) %>%
  group_by(sin) %>% mutate(N=n()) %>% ungroup %>% filter(N > 2) %>%
  select(oid, sin, date, easting, northing, all_of(hab_i$full)) %>%
  left_join(sites, by=c("sin", "date")) %>%
  mutate(east=if_else(is.na(east), easting, east),
         north=if_else(is.na(north), northing, north)) %>%
  rename(obsid=oid) %>%
  group_by(sin) %>% mutate(lon=median(east), lat=median(north)) %>% ungroup %>%
  mutate(siteid=as.numeric(factor(sin))) %>%
  rename(all_of(setNames(hab_i$full, hab_i$abbr))) %>%
  select(obsid, lon, lat, sin, siteid, date, all_of(hab_i$abbr)) %>%
  arrange(siteid, date)
site_hab.df <- fsa.df %>%  
  select(siteid, sin, lon, lat) %>%
  group_by(siteid) %>% slice_head(n=1) %>% ungroup
saveRDS(site_hab.df, "data/site_hab_df.rds")
fsa.df <- fsa.df %>% select(-lon, -lat, -sin)
saveRDS(fsa.df, "data/0_init/fsa_df.rds")

cefas.df <- fromJSON(readLines(url(urls$cefas), warn=F)) %>% as_tibble %>% 
  filter(!is.na(date_collected) & sin != "-99") %>%
  mutate(datetime_collected=as_datetime(date_collected),
         date=date(datetime_collected)) %>%
  filter(date > dateStart) %>%
  mutate(across(any_of(tox_i$full), ~if_else(.x == -99, NA_real_, .x)),
         across(any_of(tox_i$full), ~if_else(.x < 0, 0, .x))) %>%
  group_by(sin, date) %>% slice_head(n=1) %>% ungroup %>%
  group_by(sin) %>% mutate(N=n()) %>% ungroup %>% filter(N > 2) %>%
  select(oid, sin, date, easting, northing, all_of(tox_i$full)) %>%
  left_join(sites, by=c("sin", "date")) %>%
  mutate(east=if_else(is.na(east), easting, east),
         north=if_else(is.na(north), northing, north)) %>%
  rename(obsid=oid) %>%
  group_by(sin) %>% mutate(lon=median(east), lat=median(north)) %>% ungroup %>%
  filter(lat > 500000) %>%
  mutate(siteid=as.numeric(factor(sin))) %>%
  rename(all_of(setNames(tox_i$full, tox_i$abbr))) %>%
  select(obsid, lon, lat, sin, siteid, date, all_of(tox_i$abbr)) %>%
  arrange(siteid, date)
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
        "pp" # Net primary production of biomass expressed as carbon per unit volume
  ),
  source=c("Reanalysis", "Analysis&Forecast")) %>%
  mutate(server=if_else(source=="Reanalysis", "my.cmems-du.eu", "nrt.cmems-du.eu"), 
         doi=glue("https://doi.org/10.48670/moi-0005{if_else(source=='Reanalysis', 8, 6)}"),
         ID=glue("cmems_mod_nws_bgc-{var}_", 
                 "{if_else(source=='Reanalysis', 'my', 'anfc')}_7km-3D_P1D-m"))
write_csv(cmems_i, "data/cmems_i.csv")

cmems_cred <- readRDS("data/cmems_cred.rds")
fsa.df <- readRDS("data/0_init/fsa_df.rds")
cefas.df <- readRDS("data/0_init/cefas_df.rds")
get_CMEMS(userid=cmems_cred$userid, pw=cmems_cred$pw, 
          i.df=cmems_i, bbox=UK_bbox, 
          nDays_buffer=nDays_avg, dateRng=range(c(fsa.df$date, cefas.df$date)), 
          out.dir="data/0_init/cmems/")

rean.f <- dirf("data/0_init/cmems", "cmems.*Reanalysis.rds")
anfo.f <- dirf("data/0_init/cmems", "cmems.*Forecast.rds")
cmems.df <- bind_rows(
  readRDS(rean.f[1]) %>%
    bind_cols(map_dfc(rean.f[-1], ~readRDS(.x) %>% select(5))),
  readRDS(anfo.f[1]) %>%
    bind_cols(map_dfc(anfo.f[-1], ~readRDS(.x) %>% select(5)))
  ) %>%
  arrange(date, lon, lat) %>%
  group_by(date, source) %>%
  mutate(cmems_id=row_number()) %>%
  group_by(date, cmems_id, lon, lat) %>%
  summarise(across(any_of(cmems_i$var), ~mean(.x, na.rm=T))) %>%
  ungroup
  mutate(chl=log1p(chl),
         kd=log(kd),
         no3=log1p(no3),
         o2=log(o2),
         phyc=log1p(phyc),
         po4=log1p(po4)) 
saveRDS(cmems.df, glue("data/0_init/cmems_end_{max(cmems.df$date)}.rds"))



# WeStCOMS-WRF ------------------------------------------------------------

# There are three different resolutions / domains, from coarse but expansive
# (d01) to finer but restricted (d03). The getWRF() function nests all three, 
# selecting the highest resolution. Setting wrf.dir="https" will download the
# files from the public SAMS THREDDS server.
#  - Wind speed
#  - Wind direction
#  - Shortwave radiation
#  - Precipitation
#  - Sea surface temperature (not used due to extensive NAs)

fsa.df <- readRDS("data/0_init/fsa_df.rds")
cefas.df <- readRDS("data/0_init/cefas_df.rds")
wrf.dir <- ifelse(.Platform$OS.type=="unix",
                  "https",#"/media/archiver/common/sa01da-work/WRF/Archive/",
                  "D:/hydroOut/WRF/Archive/")
get_WRF(wrf.dir=wrf.dir, nDays_buffer=nDays_avg, 
        # dateRng=c(ymd("2016-01-07"), max(c(fsa.df$date, cefas.df$date))),
        dateRng=range(c(fsa.df$date, cefas.df$date)), 
        out.dir="data/0_init/")

# read and subset WRF domains to nest higher res within lower res
wrf.out <- "data/0_init/wrf/"
wrf.df <- subset_WRF("d01", wrf.out, v2_start="2019-04-01") %>%
  bind_rows(subset_WRF("d02", wrf.out, v2_start="2019-04-01")) %>%
  bind_rows(subset_WRF("d03", wrf.out, v2_start="2019-04-01")) %>%
  filter(!is.na(date)) %>%
  arrange(date, res, i) %>%
  group_by(date) %>%
  mutate(wrf_id=row_number()) %>%
  ungroup %>%
  mutate(Shortwave=log1p(Shortwave),
         Precip=log1p(pmax(Precip, 0)*3600*24*1000), # m/s to mm/day
         UV=log1p(UV))
saveRDS(wrf.df, glue("data/0_init/wrf_end_{max(wrf.df$date)}.rds"))



# pairwise distances ------------------------------------------------------

# HABs
site_hab.df <- readRDS("data/site_hab_df.rds")
path.ls <- get_shortestPaths(ocean.path="data/northAtlantic_footprint.tif", 
                             site.df=site_hab.df, 
                             transMx.path="data/mesh_tmx.rds", recalc_transMx=F)
write_csv(path.ls$dist.df, "data/site_hab_pairwise_distances.csv")
site_hab.df <- path.ls$site.df # slightly modified lat/lon to fit within ocean mesh
saveRDS(site_hab.df, "data/site_hab_df.rds")
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

# toxins
site_tox.df <- readRDS("data/site_tox_df.rds")
path.ls <- get_shortestPaths(ocean.path="data/northAtlantic_footprint.tif", 
                             site.df=site_tox.df, 
                             transMx.path="data/mesh_tmx.rds", recalc_transMx=F)
write_csv(path.ls$dist.df, "data/site_tox_pairwise_distances.csv")
site_tox.df <- path.ls$site.df # slightly modified lat/lon to fit within ocean mesh
saveRDS(site_tox.df, "data/site_tox_df.rds")
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



# fetch and bearing -------------------------------------------------------

# Wave fetch and bearing with the most open water
# https://doi.org/10.6084/m9.figshare.12029682.v1

site_hab.df <- readRDS("data/site_hab_df.rds")
site_hab.df <- site_hab.df %>%
  get_fetch(., "data/log10_eu200m1a.tif") %>%
  get_openBearing(., "data/northAtlantic_footprint.gpkg", buffer=200e3)
saveRDS(site_hab.df, "data/site_hab_df.rds")

site_tox.df <- readRDS("data/site_tox_df.rds")
site_tox.df <- site_tox.df %>%
  get_fetch(., "data/log10_eu200m1a.tif") %>%
  get_openBearing(., "data/northAtlantic_footprint.gpkg", buffer=200e3)
saveRDS(site_tox.df, "data/site_tox_df.rds")



# autoregressive terms ----------------------------------------------------

# Observed densities and bloom states for focal HAB taxa and toxins
# * = [t]
# *1 = [t-1]
# *2 = [t-2]
# N: reported density or concentration
# lnN: ln(N + 1)
# tl: HABReports traffic light for density or concentration N
# alert: HABReports action (0_none, 1_warn, 2_alert)
# lnNAvg: regional average of lnN within previous week
# prAlertAvg: regional average of alerting sites within previous week

hab.df <- calc_y_features(readRDS("data/0_init/fsa_df.rds"),
                          hab_i, tl_hab,
                          readRDS("data/site_hab_neighbors_100km.rds"))
saveRDS(hab.df, "data/0_init/hab_obs.rds")

tox.df <- calc_y_features(readRDS("data/0_init/cefas_df.rds"),
                          tox_i, tl_tox,
                          readRDS("data/site_tox_neighbors_100km.rds"))
saveRDS(tox.df, "data/0_init/tox_obs.rds")



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
st_write(site_hab.sf, "data/site_hab_sf.gpkg", append=F)

site_tox.df <- readRDS("data/site_tox_df.rds")
site_tox.sf <- site_tox.df %>% 
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  select(-sin) %>%
  st_buffer(dist=100e3) %>%
  split_to_NSEW()
st_write(site_tox.sf, "data/site_tox_sf.gpkg", append=F)



# HAB status for toxins ---------------------------------------------------

# Calculate average HAB densities surrounding each cefas site

habAvg_tox.df <- summarise_hab_states(
  site_tox.sf=st_read("data/site_tox_sf.gpkg") %>%
    group_by(siteid) %>% summarise(), 
  site_hab.sf=readRDS("data/site_hab_df.rds") %>% 
    select(siteid, lon, lat) %>% st_as_sf(coords=c("lon", "lat"), crs=27700), 
  tox.obs=readRDS("data/0_init/cefas_df.rds") %>% select(obsid, siteid, date), 
  hab.df=readRDS("data/0_init/hab_obs.rds")
)
saveRDS(habAvg_tox.df, "data/0_init/tox_habAvg.rds")



# extract sites -----------------------------------------------------------

# . CMEMS  site:date ------------------------------------------------------
cmems_i <- list(all=c("chl", "kd", "no3", "o2", "ph", "phyc", "po4", "pp"))
cmems.df <- readRDS(last(dirf("data/0_init", "cmems_end.*rds")))
cmems.sf <- cmems.df %>% select(date, lon, lat, cmems_id) %>% 
  filter(date==min(date)) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326)

# HABs
site_hab.df <- readRDS("data/site_hab_df.rds") %>% select(-cmems_id)
site_hab.df <- site_hab.df %>% find_nearest_feature_id(cmems.sf, "cmems_id")
saveRDS(site_hab.df, "data/site_hab_df.rds")
cmems.site_hab <- extract_env_pts(site_hab.df, cmems_i$all, cmems.df %>% mutate(version=1), cmems_id, "cmems_id")
saveRDS(cmems.site_hab, "data/0_init/cmems_sitePt_hab.rds")

# toxins
site_tox.df <- readRDS("data/site_tox_df.rds") %>% select(-cmems_id)
site_tox.df <- site_tox.df %>% find_nearest_feature_id(cmems.sf, "cmems_id")
saveRDS(site_tox.df, "data/site_tox_df.rds")
cmems.site_tox <- extract_env_pts(site_tox.df, cmems_i$all, cmems.df %>% mutate(version=1), cmems_id, "cmems_id")
saveRDS(cmems.site_tox, "data/0_init/cmems_sitePt_tox.rds")


# . CMEMS  buffer:date ----------------------------------------------------
cmems_i <- list(all=c("chl", "kd", "no3", "o2", "ph", "phyc", "po4", "pp"))
cmems.df <- readRDS(last(dirf("data/0_init", "cmems_end.*rds")))
cmems.sf <- cmems.df %>% select(date, lon, lat, cmems_id) %>% 
  filter(date==min(date)) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326)

# HABs
site.buffer_hab <- st_read("data/site_hab_sf.gpkg") %>%
  find_buffer_intersect_ids(., cmems.sf, "cmems_id")
cmems.buffer_hab <- extract_env_buffers(site.buffer_hab, cmems_i, cmems.df, "cmems_id")
saveRDS(cmems.buffer_hab, "data/0_init/cmems_siteBufferNSEW_hab.rds")

# toxins
site.buffer_tox <- st_read("data/site_tox_sf.gpkg") %>%
  find_buffer_intersect_ids(., cmems.sf, "cmems_id")
cmems.buffer_tox <- extract_env_buffers(site.buffer_tox, cmems_i, cmems.df, "cmems_id")
saveRDS(cmems.buffer_tox, "data/0_init/cmems_siteBufferNSEW_tox.rds")


# . WRF  site:date --------------------------------------------------------
wrf_i <- list(all=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea=c("U", "V", "UV", "Shortwave", "Precip"),
              land=c("sst"))
wrf_versions <- map(seq_along(dir("data/0_init/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/0_init/wrf", glue("domain_d0._{.x}")), readRDS) %>%
                      arrange(res, i) %>%
                      mutate(wrf_id=row_number()) %>%
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(last(dirf("data/0_init/", "wrf_end_.*rds"))) %>%
  mutate(sst=if_else(sst > -100, sst, NA_real_)) %>%
  rename(U=U_mn, V=V_mn, UV=UV_mn)

# HABs
site_hab.df <- readRDS("data/site_hab_df.rds") %>% select(-starts_with("wrf_id"))
site_hab.df <- map(wrf_versions, ~site_hab.df %>% find_nearest_feature_id(., .x, "wrf_id")) %>%
  reduce(full_join, by=names(site_hab.df), suffix=paste0(".", seq_along(wrf_versions)))
saveRDS(site_hab.df, "data/site_hab_df.rds")
site_hab.versions <- grep("wrf_id", names(site_hab.df), value=T)
wrf.site_hab <- extract_env_pts(site_hab.df, wrf_i$all, wrf.df, wrf_id, site_hab.versions)
saveRDS(wrf.site_hab, "data/0_init/wrf_sitePt_hab.rds")

# toxins
site_tox.df <- readRDS("data/site_tox_df.rds") %>% select(-starts_with("wrf_id"))
site_tox.df <- map(wrf_versions, ~site_tox.df %>% find_nearest_feature_id(., .x, "wrf_id")) %>%
  reduce(full_join, by=names(site_tox.df), suffix=paste0(".", seq_along(wrf_versions)))
saveRDS(site_tox.df, "data/site_tox_df.rds")
site_tox.versions <- grep("wrf_id", names(site_tox.df), value=T)
wrf.site_tox <- extract_env_pts(site_tox.df, wrf_i$all, wrf.df, wrf_id, site_tox.versions)
saveRDS(wrf.site_tox, "data/0_init/wrf_sitePt_tox.rds")


# . WRF  buffer:date ------------------------------------------------------
wrf_i <- list(all=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea=c("U", "V", "UV", "Shortwave", "Precip"),
              land=c("sst"))
wrf_versions <- map(seq_along(dir("data/0_init/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/0_init/wrf", glue("domain_d0._{.x}")), readRDS) %>%
                      arrange(res, i) %>%
                      mutate(wrf_id=row_number()) %>%
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(last(dirf("data/0_init/", "wrf_end_.*rds"))) %>%
  mutate(sst=if_else(sst > -100, sst, NA_real_)) %>%
  rename(U=U_mn, V=V_mn, UV=UV_mn)

# HABs
site.buffer_hab <- map(wrf_versions, 
                       ~st_read("data/site_hab_sf.gpkg") %>% 
                         find_buffer_intersect_ids(., .x, "wrf_id")) %>%
  reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer_hab <- extract_env_buffers(site.buffer_hab, wrf_i, wrf.df, paste0("wrf_id.", 1:2))
saveRDS(wrf.buffer_hab, "data/0_init/wrf_siteBufferNSEW_hab.rds")

# toxins
site.buffer_tox <- map(wrf_versions, 
                       ~st_read("data/site_tox_sf.gpkg") %>% 
                         find_buffer_intersect_ids(., .x, "wrf_id")) %>%
  reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer_tox <- extract_env_buffers(site.buffer_tox, wrf_i, wrf.df, paste0("wrf_id.", 1:2))
saveRDS(wrf.buffer_tox, "data/0_init/wrf_siteBufferNSEW_tox.rds")



# compile -----------------------------------------------------------------

# HABs
hab.ls <- load_datasets("0_init", "hab")
hab.df <- hab.ls$site %>% select(-sin) %>%
  right_join(hab.ls$obs, by="siteid", multiple="all") %>%
  left_join(hab.ls$cmems.pt, by=c("cmems_id", "date")) %>%
  left_join(hab.ls$cmems.buf %>%
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) %>%
  select(-wrf_id.1, -wrf_id.2, -version) %>%
  left_join(hab.ls$wrf.pt %>% select(-version), by=c("wrf_id", "date")) %>%
  left_join(hab.ls$wrf.buf %>%
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  mutate(year=year(date),
         yday=yday(date)) %>%
  select(-contains("sst")) 
saveRDS(hab.df, "data/0_init/data_hab_all_withNA.rds")

# toxins
tox.ls <- load_datasets("0_init", "tox")
tox.df <- tox.ls$site %>% select(-sin) %>%
  right_join(tox.ls$obs, by="siteid", multiple="all") %>%
  left_join(tox.ls$habAvg %>% select(-date, -siteid), by="obsid") %>%
  left_join(tox.ls$cmems.pt, by=c("cmems_id", "date")) %>%
  left_join(tox.ls$cmems.buf %>%
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) %>%
  select(-wrf_id.1, -wrf_id.2, -version) %>%
  left_join(tox.ls$wrf.pt %>% select(-version), by=c("wrf_id", "date")) %>%
  left_join(tox.ls$wrf.buf %>%
              pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
            by=c("siteid", "date")) %>%
  mutate(year=year(date),
         yday=yday(date)) %>%
  select(-contains("sst")) 
saveRDS(tox.df, "data/0_init/data_tox_all_withNA.rds")

# variable names
grep("cmems_id|date|siteid|version", 
     unique(c(names(hab.ls$cmems.pt), 
              names(hab.ls$cmems.buf %>%
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")))),
     value=T, invert=T) %>% 
  sort %>%
  saveRDS("data/cmems_vars.rds")
grep("wrf_id|date|siteid|version", 
     unique(c(names(hab.ls$wrf.pt), 
              names(hab.ls$wrf.buf %>%
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")))),
     value=T, invert=T) %>% 
  sort %>%
  saveRDS("data/wrf_vars.rds")


