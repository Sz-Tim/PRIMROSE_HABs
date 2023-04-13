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




# sampling locations and dates --------------------------------------------

fsa.df <- fromJSON("data/0_init/copy_fsa.txt") %>% 
  as_tibble %>% 
  filter(!is.na(date_collected)) %>%
  filter(karenia_mikimotoi >= 0) %>% # -99 in karenia...?
  mutate(datetime_collected=as_datetime(date_collected),
         date=date(datetime_collected)) %>%
  filter(datetime_collected >= "2013-07-20") %>% 
  mutate(siteid=as.numeric(factor(sin))) %>%
  rename(obsid=oid) %>%
  group_by(sin) %>% 
  mutate(lon=median(easting), lat=median(northing)) %>%
  ungroup %>%
  select(-geom, -easting, -northing, -tide, -datetime_collected, -samplemethod, 
         -depth, -site, -area)

site.df <- fsa.df %>%  
  select(siteid, sin, lon, lat) %>%
  group_by(siteid) %>% slice_head(n=1) %>% ungroup
site.4326 <- site.df %>% 
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  st_buffer(dist=50e3) %>%
  st_transform(4326)




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
write_csv(cmems_i, "data/0_init/cmems_i.csv")

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
  filter(complete.cases(.))
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
get_shortestPaths(ocean.path="data/0_init/northAtlantic_footprint.tif", 
                  site.df=site.df, 
                  transMx.path="data/0_init/mesh_tmx.rds", recalc_transMx=F) %>%
  write_csv("data/0_init/fsa_site_pairwise_distances.csv")




# Fetch and bearing -------------------------------------------------------

# Wave fetch and bearing with the most open water
# https://doi.org/10.6084/m9.figshare.12029682.v1
site.df <- site.df %>%
  get_fetch(., "https://figshare.com/ndownloader/files/22107477") %>%
  get_openBearing(., "data/0_init/northAtlantic_footprint.gpkg", buffer=200e3)




# HAB densities -----------------------------------------------------------

# Observed densities and bloom states for focal taxa
# N[t]
# N[t-1]
# N[t-2]






# Compile and save --------------------------------------------------------

# Partition datasets for testing/training
