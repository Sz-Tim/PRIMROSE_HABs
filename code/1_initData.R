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
  mutate(N_ord=LETTERS[row_number()],
         alert=case_when(is.na(alert)~0,
                         alert=="warn"~1,
                         alert=="alert"~2)) %>%
  select(sp, min_ge, alert, tl, N_ord)
hab.tl <- hab.tl %>% 
  bind_rows(hab.tl %>% filter(sp=="pseudo_nitzschia_sp") %>% 
              mutate(sp="pseudo_nitzschia_delicatissima_group")) %>% 
  bind_rows(hab.tl %>% filter(sp=="pseudo_nitzschia_sp") %>% 
              mutate(sp="pseudo_nitzschia_seriata_group")) %>%
  mutate(sp=sp_i$abbr[match(sp, sp_i$full)])



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
saveRDS(site.df, "data/site_df.rds")

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

get_CMEMS(userid="tszewczyk", pw="xaj5kba*haq_TYD4pzb", 
          i.df=cmems_i, bbox=UK_bbox, 
          nDays_buffer=nDays_avg, dateRng=range(fsa.df$date), 
          out.dir="data/0_init/cmems/")

rean.f <- dir("data/0_init/cmems", "cmems.*Reanalysis.rds", full.names=T)
anfo.f <- dir("data/0_init/cmems", "cmems.*Forecast.rds", full.names=T)
cmems.df <- bind_rows(
  readRDS(rean.f[1]) %>%
    bind_cols(map_dfc(rean.f[-1], ~readRDS(.x) %>% select(5))) %>%
    filter(date < "2019-05-01"),
  readRDS(anfo.f[1]) %>%
    bind_cols(map_dfc(anfo.f[-1], ~readRDS(.x) %>% select(5)))
  ) %>%
  filter(complete.cases(.)) %>%
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
wrf.f <- dir("data/0_init/wrf/", "wrf.*rds", full.names=T)
wrf.df <- map_dfr(wrf.f, readRDS)
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
saveRDS(site.df, "data/site_df.rds")





# site buffers ------------------------------------------------------------

# Buffers for averaging environmental conditions
# CMEMS is coarse and does not resolve lochs, so some sites are not covered
# The buffer size was determined by trial and error to reduce averaging while
# ensuring all sites are represented.

coast.sf <- st_read("data/northAtlantic_footprint.gpkg")
site.sf <- map(c(0.5, 1, 5, 10, 25, 50), 
               ~site.df %>% 
                 st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
                 select(-sin) %>% mutate(buffer_km=.x) %>%
                 st_buffer(dist=.x*1e3))




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

hab.df <- fsa.df %>% 
  pivot_longer(any_of(sp_i$abbr), names_to="sp", values_to="N") %>%
  filter(!is.na(N)) %>%
  mutate(lnN=log1p(N)) %>%
  get_trafficLights(N, hab.tl) %>%
  arrange(sp, siteid, date) %>%
  group_by(sp, siteid) %>%
  get_lags(lnN, alert, date, n=2) %>%
  mutate(lnDayLag1=log(as.numeric(date-date1)), 
         lnDayLag2=log(as.numeric(date-date2)))
for(i in 1:nrow(hab.df)) {
  site_i <- hab.df$siteid[i]
  date_i <- hab.df$date[i]
  sp_i <- hab.df$sp[i]
  wk.df <- hab.df %>% select(siteid, date, lnN1, lnN2, alert1, alert2) %>%
    filter(siteid %in% site.100k$dest_c[site.100k$origin==site_i][[1]]) %>%
    filter(date <= date_i & date > date_i-7) 
  hab.df$lnNAvg1[i] <- mean(wk.df$lnN1, na.rm=T)
  hab.df$lnNAvg2[i] <- mean(wk.df$lnN2, na.rm=T)
  hab.df$prAlertAvg1[i] <- mean(wk.df$alert1 != "0_none", na.rm=T)
  hab.df$prAlertAvg2[i] <- mean(wk.df$alert2 != "0_none", na.rm=T)
  if(i %% 100 == 0) {cat(i, "of", nrow(hab.df), "\n")}
}
saveRDS(hab.df, "data/0_init/hab_densities.rds")





# Compile and save --------------------------------------------------------

# load, extract, and compile
# - Algal densities
# - CMEMS
# - WRF
# - fetch & bearing
# - calculate cos(dir)'s

# Extract CMEMS data for each site:date
cmems.df <- readRDS(dir("data/0_init", "cmems_.*rds", full.names=T))
cmems.sf <- cmems.df %>% select(date, lon, lat, cmems_id) %>% 
  filter(date==min(date)) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326)
site.df <- site.df %>%
  st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) %>%
  st_transform(4326) %>%
  mutate(cmems_id=st_nearest_feature(., cmems.sf)) %>%
  st_drop_geometry()
cmems.df <- cmems.df %>% 
  filter(cmems_id %in% site.df$cmems_id) %>%
  group_by(cmems_id) %>%
  mutate(across(any_of(unique(cmems_i$var)), 
                ~zoo::rollmean(.x, k=7, na.pad=T),
                .names="{.col}Wk")) %>%
  ungroup %>%
  mutate(day_cumul=as.numeric(date - ymd("2010-01-01")),
         yday=yday(date)) %>%
  group_by(cmems_id) %>%
  mutate(across(any_of(unique(cmems_i$var)),
                ~detrend_loess(yday, .x, span=0.2), 
                .names="{.col}Dt")) %>%
  ungroup

# Extract WRF data for each site:date


obs.df <- readRDS("data/0_init/hab_densities.rds") %>%
  left_join(readRDS("data/site_df.rds"))

# Partition datasets for testing/training
