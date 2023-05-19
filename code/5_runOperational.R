# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Operational runs



# TODO: Make a lot of functions or subscripts that can be called as needed
# Need to download new datasets, process them, and merge them without overwriting
# anything that shouldn't be overwritten, This will take some thought to avoid
# lots of repetitive code and clunky structure. 
# Refitting models can be more manual
# All new data just needs to be to be for a testing dataset

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

daysBuffer <- 10
daysForecast <- 7
UK_bbox <- list(xmin=-11, xmax=3, ymin=49, ymax=61.5)
urls <- c(fsa="fsa_counts",
          fsa_sites="fsa_sites", 
          cefas="cefas_counts",
          cefas_sites="cefas_sites") %>%
  map(~glue("http://varro:3001/{.x}"))

# NB: toxin thresholds are adjusted for ASP, AZAs, YTXs because of the extremely 
# limited occurrence of TL2 and TL3. Models thus predict only presence/absence
hab_i <- read_csv("data/i_hab.csv")
tox_i <- read_csv("data/i_tox.csv")

old_end <- readRDS("data/1_current/obs_end.rds") %>%
  map(~ymd(.x)-daysBuffer)

# sampling locations and dates --------------------------------------------

site_hab.df <- readRDS("data/site_hab_df.rds")
fsa_1 <- readRDS("data/1_current/fsa_df.rds") %>% 
  filter(date < old_end$hab, !is.na(obsid))
fsa_2 <- read_and_clean_fsa(urls$fsa, hab_i, fsa_sites, old_end$hab) %>%
  left_join(site_hab.df %>% select(sin, siteid))
if(any(is.na(fsa_2$siteid))) {
  new_sin <- paste(unique(filter(fsa_2, is.na(siteid))$sin), collapse="\n")
  stop(c("New HAB sin found\n", new_sin))
}
fsa.df <- bind_rows(fsa_1, fsa_2 %>% select(-lon, -lat)) %>%
  bind_rows(site_hab.df %>% select(siteid, sin) %>% mutate(date=max(fsa_2$date)+daysForecast))
saveRDS(fsa.df, "data/2_new/fsa_df.rds")

site_tox.df <- readRDS("data/site_tox_df.rds")
cefas_1 <- readRDS("data/1_current/cefas_df.rds") %>% 
  filter(date < old_end$tox, !is.na(obsid))
cefas_2 <- read_and_clean_cefas(urls$cefas, tox_i, cefas_sites, old_end$tox) %>%
  left_join(site_tox.df %>% select(sin, siteid))
if(any(is.na(cefas_2$siteid))) {
  new_sin <- paste(unique(filter(cefas_2, is.na(siteid))$sin), collapse="\n")
  stop(c("New toin sin found\n", new_sin))
}
cefas.df <- bind_rows(cefas_1, cefas_2 %>% select(-lon, -lat)) %>%
  bind_rows(site_tox.df %>% select(siteid, sin) %>% mutate(date=max(cefas_2$date)+daysForecast))
saveRDS(cefas.df, "data/2_new/cefas_df.rds")



# refresh datasets --------------------------------------------------------

# . CMEMS  download -------------------------------------------------------
cmems_cred <- readRDS("data/cmems_cred.rds")
cmems_i <- read_csv("data/cmems_i.csv") %>% filter(source=="Analysis&Forecast")
get_CMEMS(cmems_cred$userid, cmems_cred$pw, cmems_i, UK_bbox, daysBuffer, 
          c(old_end$cmems, Sys.Date()+7), "data/00_env/cmems/")
anfo.f <- dirf("data/00_env/cmems", "cmems.*Forecast.rds")
cmems.df <- map(anfo.f, readRDS) %>% 
  reduce(full_join) %>%
  arrange(date, lon, lat) %>%
  na.omit() %>%
  group_by(date) %>%
  mutate(cmems_id=row_number()) %>%
  ungroup %>%
  mutate(chl=log1p(chl),
         kd=log(kd),
         no3=log1p(no3),
         o2=log(o2),
         phyc=log1p(phyc),
         po4=log1p(po4)) 
saveRDS(cmems.df, glue("data/2_new/cmems_end_{max(cmems.df$date)}.rds"))

# . CMEMS  site:date ------------------------------------------------------
# replace previous forecasts with recent observations
# add new forecasts for predictions
cmems_i <- list(all=c("chl", "kd", "no3", "o2", "ph", "phyc", "po4", "pp"))
cmems.df <- readRDS(last(dirf("data/2_new", "cmems_end.*rds"))) %>% mutate(version=1)
cmems.sf <- cmems.df %>% 
  filter(date==min(date)) %>%
  select(date, lon, lat, cmems_id) %>%
  st_as_sf(., coords=c("lon", "lat"), crs=4326)

# HABs
site_hab.df <- readRDS("data/site_hab_df.rds") 
cmems.site_hab_1 <- readRDS("data/1_current/cmems_sitePt_hab.rds")
cmems.site_hab_2 <- extract_env_pts(site_hab.df, cmems_i$all, cmems.df, cmems_id, "cmems_id") %>%
  na.omit
cmems.site_hab <- bind_rows(cmems.site_hab_1 %>% filter(date < min(cmems.site_hab_2$date)),
                            cmems.site_hab_2)
saveRDS(cmems.site_hab %>% na.omit, "data/2_new/cmems_sitePt_hab.rds")

# toxins
site_tox.df <- readRDS("data/site_tox_df.rds") 
cmems.site_tox_1 <- readRDS("data/1_current/cmems_sitePt_tox.rds")
cmems.site_tox_2 <- extract_env_pts(site_tox.df, cmems_i$all, cmems.df, cmems_id, "cmems_id") %>%
  na.omit
cmems.site_tox <- bind_rows(cmems.site_tox_1 %>% filter(date < min(cmems.site_tox_2$date)),
                            cmems.site_tox_2)
saveRDS(cmems.site_tox %>% na.omit, "data/2_new/cmems_sitePt_tox.rds")


# . CMEMS  buffer:date ----------------------------------------------------
# HABs
site.buffer_hab <- st_read("data/site_hab_sf.gpkg") %>%
  find_buffer_intersect_ids(., cmems.sf, "cmems_id")
cmems.buffer_hab_1 <- readRDS("data/1_current/cmems_siteBufferNSEW_hab.rds")
cmems.buffer_hab_2 <- extract_env_buffers(site.buffer_hab, cmems_i, cmems.df, "cmems_id") %>%
  na.omit
cmems.buffer_hab <- bind_rows(cmems.buffer_hab_1 %>% filter(date < min(cmems.buffer_hab_2$date)),
                              cmems.buffer_hab_2)
saveRDS(cmems.buffer_hab, "data/2_new/cmems_siteBufferNSEW_hab.rds")

# toxins
site.buffer_tox <- st_read("data/site_tox_sf.gpkg") %>%
  find_buffer_intersect_ids(., cmems.sf, "cmems_id")
cmems.buffer_tox_1 <- readRDS("data/1_current/cmems_siteBufferNSEW_tox.rds")
cmems.buffer_tox_2 <- extract_env_buffers(site.buffer_tox, cmems_i, cmems.df, "cmems_id")
cmems.buffer_tox <- bind_rows(cmems.buffer_tox_1 %>% filter(date < min(cmems.buffer_tox_2$date)),
                              cmems.buffer_tox_2)
saveRDS(cmems.buffer_tox, "data/2_new/cmems_siteBufferNSEW_tox.rds")


# . WRF  download ---------------------------------------------------------
wrf.out <- "data/00_env/wrf/"
get_WRF(wrf.dir="https", nDays_buffer=daysBuffer, 
        dateRng=c(old_end$wrf, Sys.Date()+7), out.dir=wrf.out)
get_WRF(wrf.dir="https", nDays_buffer=0, 
        dateRng=c(old_end$wrf, Sys.Date()+7), out.dir=wrf.out, 
        forecast=T)
new_obs <- dir(wrf.out, "^wrfF?_") %>% str_split_fixed("_", 2) %>% as_tibble %>% 
  group_by(V2) %>% summarise(N=n()) %>% ungroup %>% filter(N>2)
if(nrow(new_obs) > 0) {
  file.remove(glue("{wrf.out}wrfF_{new_obs$V2}"))
}
wrf.df <- aggregate_WRF(wrf.out)
saveRDS(wrf.df, glue("data/0_init/wrf_end_{max(wrf.df$date)}.rds"))

# . WRF  site:date --------------------------------------------------------
wrf_i <- list(all=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea=c("U", "V", "UV", "Shortwave", "Precip"),
              land=c("sst"))
wrf_versions <- map(seq_along(dir("data/00_env/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/00_env/wrf", glue("domain_d0._{.x}")), readRDS) %>%
                      arrange(res, i) %>%
                      mutate(wrf_id=row_number()) %>%
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(last(dirf("data/2_new/", "wrf_end_.*rds"))) 

# HABs
site_hab.df <- readRDS("data/site_hab_df.rds")
site_hab.versions <- grep("wrf_id", names(site_hab.df), value=T)
wrf.site_hab <- extract_env_pts(site_hab.df, wrf_i$all, wrf.df, wrf_id, site_hab.versions)
saveRDS(wrf.site_hab, "data/2_new/wrf_sitePt_hab.rds")

# toxins
site_tox.df <- readRDS("data/site_tox_df.rds")
site_tox.versions <- grep("wrf_id", names(site_tox.df), value=T)
wrf.site_tox <- extract_env_pts(site_tox.df, wrf_i$all, wrf.df, wrf_id, site_tox.versions)
saveRDS(wrf.site_tox, "data/2_new/wrf_sitePt_tox.rds")


# . WRF  buffer:date ------------------------------------------------------
wrf_i <- list(all=c("U", "V", "UV", "Shortwave", "Precip", "sst"),
              sea=c("U", "V", "UV", "Shortwave", "Precip"),
              land=c("sst"))
wrf_versions <- map(seq_along(dir("data/00_env/wrf", "^domain_d01")), 
                    ~map_dfr(dirf("data/00_env/wrf", glue("domain_d0._{.x}")), readRDS) %>%
                      arrange(res, i) %>%
                      mutate(wrf_id=row_number()) %>%
                      st_as_sf(coords=c("lon", "lat"), remove=F, crs=4326))
wrf.df <- readRDS(last(dirf("data/2_new/", "wrf_end_.*rds"))) 

# HABs
site.buffer_hab <- map(wrf_versions, 
                       ~st_read("data/site_hab_sf.gpkg") %>% 
                         find_buffer_intersect_ids(., .x, "wrf_id")) %>%
  reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer_hab <- extract_env_buffers(site.buffer_hab, wrf_i, wrf.df, paste0("wrf_id.", 1:2))
saveRDS(wrf.buffer_hab, "data/2_new/wrf_siteBufferNSEW_hab.rds")

# toxins
site.buffer_tox <- map(wrf_versions, 
                       ~st_read("data/site_tox_sf.gpkg") %>% 
                         find_buffer_intersect_ids(., .x, "wrf_id")) %>%
  reduce(full_join, by=c("siteid", "quadrant"), suffix=paste0(".", seq_along(wrf_versions)))
wrf.buffer_tox <- extract_env_buffers(site.buffer_tox, wrf_i, wrf.df, paste0("wrf_id.", 1:2))
saveRDS(wrf.buffer_tox, "data/2_new/wrf_siteBufferNSEW_tox.rds")

# calculate y
hab.df <- calc_y_features(readRDS("data/2_new/fsa_df.rds"),
                          hab_i, tl_hab,
                          readRDS("data/site_hab_neighbors_100km.rds"))
saveRDS(hab.df, "data/2_new/hab_obs.rds")

tox.df <- calc_y_features(readRDS("data/2_new/cefas_df.rds"),
                          tox_i, tl_tox,
                          readRDS("data/site_tox_neighbors_100km.rds"))
saveRDS(tox.df, "data/2_new/tox_obs.rds")

# hab avg
habAvg_tox.df <- summarise_hab_states(
  site_tox.sf=st_read("data/site_tox_sf.gpkg") %>%
    group_by(siteid) %>% summarise(), 
  site_hab.sf=readRDS("data/site_hab_df.rds") %>% 
    select(siteid, lon, lat) %>% st_as_sf(coords=c("lon", "lat"), crs=27700), 
  tox.obs=readRDS("data/2_new/cefas_df.rds") %>% select(obsid, siteid, date), 
  hab.df=readRDS("data/2_new/hab_obs.rds")
)
saveRDS(habAvg_tox.df, "data/2_new/tox_habAvg.rds")

# compile -----------------------------------------------------------------

# HABs
hab.ls <- load_datasets("2_new", "hab")
saveRDS(hab.ls$compiled, "data/2_new/data_hab_all.rds")

# toxins
tox.ls <- load_datasets("2_new", "tox")
saveRDS(tox.ls$compiled, "data/2_new/data_tox_all.rds")

obs_end <- list(
  hab=max(hab.ls$obs$date),
  tox=max(tox.ls$obs$date),
  cmems=ymd(str_sub(dir("data/2_new", "cmems_end"), 11, 20)),
  wrf=ymd(str_sub(dir("data/2_new", "wrf_end"), 9, 18))
)
saveRDS(obs_end, "data/2_new/obs_end.rds")

write_to_current <- F
if(write_to_current) {
  file.copy(dirf("data/2_new/", "fsa_df"), "data/1_current/")
  file.copy(dirf("data/2_new/", "cefas_df"), "data/1_current/")
  file.copy(dirf("data/2_new/", "_obs.rds"), "data/1_current/")
  file.copy(dirf("data/2_new/", "_habAvg.rds"), "data/1_current/")
  file.copy(dirf("data/2_new/", "_sitePt_"), "data/1_current/")
  file.copy(dirf("data/2_new/", "_siteBufferNSEW_"), "data/1_current/")
  file.copy(dirf("data/2_new/", "data_.*_all.rds"), "data/1_current/")
  file.copy(dirf("data/2_new/", "obs_end"), "data/1_current/")
}




# load models -------------------------------------------------------------
rm(list=ls()); gc()
pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", "nnet", "randomForest", "glmnet", 
          "brms", "bayesian", "doParallel", "foreach")
lapply(pkgs, library, character.only=T)
source("code/00_fn.R")
fit.dir <- glue("out/model_fits/{covSet}/")
cv.dir <- glue("{fit.dir}/cv/")
out.dir <- glue("out/{covSet}/")
y_i <- bind_rows(read_csv("data/i_hab.csv") %>% arrange(abbr) %>% mutate(type="hab"),
                 read_csv("data/i_tox.csv") %>% arrange(abbr) %>% mutate(type="tox"))

col_metadata <- c("obsid", "y", "date", "year", "yday", "siteid", "lon", "lat")
col_resp <- c("lnN", "tl", "alert")
col_cmems <- readRDS("data/cmems_vars.rds")
col_wrf <- readRDS("data/wrf_vars.rds")

all_covs <- list(
  spacetime=c("ydayCos", "ydaySin", "ydaySinXydayCos", 
              "latz", "lonz", "lonzXlatz"),
  main=c(
    "fetch", 
    "lnNWt1", "lnNAvg1", "prAlertAvg1", "alert1A1", 
    "lnNWt2", "lnNAvg2", "prAlertAvg2", "alert2A1",
    col_cmems, col_wrf
  ),
  interact=c(
    paste("UWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
    paste("VWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
    paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X"),
    paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X")
  ),
  hab=c(outer(filter(y_i, type=="hab")$abbr, c("lnNAvg", "prA"), "paste0"))
)

covs_exclude <- switch(covSet,
                       full="NA",
                       test="Xfetch|Dt|Delta|Avg",
                       noDt="Dt",
                       noDtDelta="Dt|Delta",
                       noInt="Xfetch",
                       noDtDeltaInt="Xfetch|Dt|Delta")

obs.ls <- map_dfr(dirf("data/1_current", "data_.*_all.rds"), readRDS) %>%
  filter(date >= (Sys.Date()-daysBuffer)) %>%
  select(all_of(col_metadata), all_of(col_resp), 
         "alert1", "alert2", any_of(unname(unlist(all_covs)))) %>%
  mutate(across(starts_with("alert"), ~factor(.x)),
         across(starts_with("tl"), ~factor(.x, ordered=T))) %>%
  group_by(y) %>%
  group_split() %>%
  map(~.x %>% select(where(~any(!is.na(.x)))) %>% na.omit)



# generate predictions ----------------------------------------------------
  for(i in seq_along(obs.ls)) {
    y <- obs.ls[[i]]$y[1]
    y_i.i <- filter(y_i, abbr==y)
    responses <- c(alert="alert", tl="tl")
    wt.ls <- readRDS(glue("{out.dir}/{y}_wt_ls.rds"))
    
    oos.ls <- map(responses, ~summarise_predictions(obs.ls[[i]], .x, fit.dir, y_i.i))
    oos.ls$alert <- oos.ls %>%
      map(~.x %>% select(y, obsid, siteid, date, ends_with("_A1"))) %>%
      reduce(full_join)
    
    oos.ls <- map(responses, ~calc_ensemble(oos.ls, wt.ls, .x, y_i.i))
    oos.ls$alert <- full_join(
      oos.ls$alert, 
      oos.ls[-which(responses=="alert")] %>%
        map(~.x %>% select(y, obsid, siteid, date, matches("ens_.*_A1"))) %>%
        reduce(full_join)
    )
    saveRDS(oos.ls, glue("{out.dir}/{y}_oos_ls.rds"))
  }


# store output ------------------------------------------------------------

