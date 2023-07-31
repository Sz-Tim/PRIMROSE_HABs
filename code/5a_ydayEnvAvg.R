# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Operational runs: Recalculate yday environmental averages



# setup -------------------------------------------------------------------
library(tidyverse); library(glue)
source("code/00_fn.R")

pt_data <- c("sitePt_hab", "sitePt_tox")
buf_data <- c("siteBufferNSEW_hab", "siteBufferNSEW_tox")



# CMEMS -------------------------------------------------------------------

walk(pt_data, 
    ~readRDS(glue("data/1_current/cmems_{.x}.rds")) |> 
      mutate(yday=yday(date)) |>
      group_by(cmems_id, version, yday) |>
      summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) |>
      ungroup() |>
      saveRDS(glue("data/1_current/ydayAvg_cmems_{.x}.rds")))
walk(buf_data, 
     ~readRDS(glue("data/1_current/cmems_{.x}.rds")) |> 
       mutate(yday=yday(date)) |>
       group_by(siteid, quadrant, yday) |>
       summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) |>
       ungroup() |>
       saveRDS(glue("data/1_current/ydayAvg_cmems_{.x}.rds")))



# WRF ---------------------------------------------------------------------

walk(pt_data, 
     ~readRDS(glue("data/1_current/wrf_{.x}.rds")) |> 
       mutate(yday=yday(date)) |>
       group_by(wrf_id, version, yday) |>
       summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) |>
       ungroup() |>
       saveRDS(glue("data/1_current/ydayAvg_wrf_{.x}.rds")))
walk(buf_data, 
     ~readRDS(glue("data/1_current/wrf_{.x}.rds")) |> 
       mutate(yday=yday(date)) |>
       group_by(siteid, quadrant, yday) |>
       summarise(across(where(is.numeric), ~mean(.x, na.rm=T))) |>
       ungroup() |>
       saveRDS(glue("data/1_current/ydayAvg_wrf_{.x}.rds")))

