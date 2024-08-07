# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Functions: data preparation


read_and_clean_sites <- function(url_sites, dateStart) {
  paste0(url_sites, "?fromdate=gte.", dateStart) |>
    url() |>
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |>
    filter(east < 7e5,
           north < 125e4,
           !(east==0 & north==0),
           sin != "-99",
           !is.na(fromdate) & !is.na(todate),
           fromdate != todate) |>
    mutate(fromdate=date(fromdate), todate=date(todate)) |>
    rowwise() |>
    mutate(date=list(seq(fromdate, todate, by=1))) |>
    ungroup() |>
    arrange(sin, fromdate) |>
    select(sin, east, north, date) |>
    unnest(date) |>
    arrange(sin, date) |>
    group_by(sin, date) |>
    slice_head(n=1) |>
    ungroup()
}



read_and_clean_fsa <- function(url_fsa, hab_i, sites, dateStart="2016-01-01") {
  paste0(url_fsa, "?date_collected=gte.", dateStart) |>
    url() |> 
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |> 
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) |>
    mutate(across(any_of(hab_i$full), ~na_if(.x, -99))) |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, all_of(hab_i$full)) |>
    left_join(sites, by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    rename(all_of(setNames(hab_i$full, hab_i$abbr))) |>
    select(obsid, lon, lat, sin, date, all_of(hab_i$abbr)) |>
    arrange(sin, date) 
}


read_and_clean_cefas <- function(url_cefas, tox_i, sites, dateStart="2016-01-01") {
  paste0(url_cefas, "?date_collected=gte.", dateStart) |>
    url() |> 
    readLines(warn=F) |>
    fromJSON() |> as_tibble() |> 
    filter(!is.na(date_collected) & sin != "-99") |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) |>
    mutate(across(any_of(tox_i$full), ~if_else(.x == -99, NA_real_, .x)),
           across(any_of(tox_i$full), ~if_else(.x < 0, 0, .x))) |>
    group_by(sin, date) |> slice_head(n=1) |> ungroup() |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, all_of(tox_i$full)) |>
    left_join(sites, by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    filter(lat > 500000) |>
    rename(all_of(setNames(tox_i$full, tox_i$abbr))) |>
    select(obsid, lon, lat, sin, date, all_of(tox_i$abbr)) |>
    arrange(sin, date)
}


read_and_clean_fish <- function(url_mowi, url_ssf, fish_i, sites, dateStart="2016-01-01") {
  bind_rows(url(paste0(url_mowi, "?date_collected=gte.", dateStart)) |> 
              readLines(warn=F) |>
              fromJSON() |> as_tibble(),
            url(paste0(url_ssf, "?date_collected=gte.", dateStart)) |> 
              readLines(warn=F) |>
              fromJSON() |> as_tibble()) |> 
    filter(!is.na(date_collected)) |>
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) |>
    mutate(across(any_of(fish_i$full), ~na_if(.x, -99))) |>
    group_by(sin) |> mutate(N=n()) |> ungroup() |> filter(N > 2) |>
    select(oid, sin, date, easting, northing, any_of(fish_i$full)) |>
    left_join(sites, by=c("sin", "date")) |>
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) |>
    rename(obsid=oid) |>
    group_by(sin) |> mutate(lon=median(east), lat=median(north)) |> ungroup() |>
    rename(any_of(setNames(fish_i$full, fish_i$abbr))) |>
    select(obsid, lon, lat, sin, date, any_of(fish_i$abbr)) |>
    arrange(sin, date) |>
    filter(sin!=0) |>
    group_by(date, sin) |> 
    summarise(across(where(is.numeric), ~mean(.x, na.rm=T)))
}





#' Download CMEMS layers
#'
#' @param userid 
#' @param pw 
#' @param i.df 
#' @param bbox 
#' @param nDays_buffer 
#' @param dateRng 
#' @param out.dir 
#'
#' @return
#' @export
#'
#' @examples
get_CMEMS <- function(userid, pw, i.df, bbox, nDays_buffer, dateRng, out.dir,
                      toolbox=TRUE) {
  if(toolbox) {
    save(list=ls(all.names=TRUE), file="temp/get_CMEMS.RData")
    system2("bash", paste0(getwd(), "/code/00_getCMEMS.sh"))
    file.remove("temp/get_CMEMS.RData")
    return("Finished running /code/getCMEMS.sh")
  } 
  
  library(tidyverse); library(ncdf4); library(lubridate); library(glue)
  
  for(i in 1:nrow(i.df)) {
    # load .nc file
    url <- paste0("https://", userid, ":", pw, "@", 
                  i.df$server[i], "/thredds/dodsC/", i.df$ID[i])
    nc <- nc_open(url)
    
    # get metadata and identify lon/lat/time indexes to extract
    nc.lon <- ncvar_get(nc, "longitude")
    nc.lat <- ncvar_get(nc, "latitude")
    nc.time <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")
    nc.origin <- as_datetime(str_split_fixed(time_units$value, " ", 3)[3])
    nc.dt <- switch(str_split_fixed(time_units$value, " ", 3)[1], 
                    seconds=dseconds,
                    minutes=dminutes,
                    hours=dhours)
    nc.date <- date(nc.origin + nc.dt(nc.time))
    lon_i <- which(between(nc.lon, bbox$xmin, bbox$xmax))
    lat_i <- which(between(nc.lat, bbox$ymin, bbox$ymax))
    time_i <- which(between(nc.date, dateRng[1]-nDays_buffer, dateRng[2]+nDays_buffer))
    var.start <- c(lon_i[1], lat_i[1], 1, time_i[1])
    var.count <- c(length(lon_i), length(lat_i), 1, length(time_i)) 
    
    # extract variable
    if(length(names(nc$var)) > 1) { cat("Too many vars in nc:", names(nc$var))}
    nc.var <- ncvar_get(nc, names(nc$var), start=var.start, count=var.count)
    nc_close(nc)
    
    # reshape and save
    expand_grid(date=nc.date[time_i],
                lat=c(nc.lat[lat_i]),
                lon=c(nc.lon[lon_i])) |>
      mutate(source=i.df$source[i],
             vals=c(nc.var)) |>
      rename_with(~gsub("vals", i.df$var[i], .x)) |>
      saveRDS(glue("{out.dir}/cmems_{i.df$var[i]}_{i.df$source[i]}.rds"))
  }
}













#' Download WRF data
#'
#' @param wrf.dir 
#' @param nDays_buffer 
#' @param dateRng 
#' @param out.dir 
#'
#' @return
#' @export
#'
#' @examples
get_WRF <- function(wrf.dir, nDays_buffer, dateRng, out.dir, forecast=F) {
  library(tidyverse); library(ncdf4); library(lubridate); library(glue); 
  library(xml2); library(rvest)
  dir.create(out.dir, showWarnings=F)
  
  # metadata for all WRF files within timespan
  if(grepl("https", wrf.dir)) {
    thredd_base <- "https://thredds.sams.ac.uk/thredds/"
    Archive <- ifelse(forecast, "Archive_forecast", "Archive")
    catalog <- ifelse(forecast, "F/catalog.html", "/catalog.html")
    
    thredds_head <- glue("{thredd_base}catalog/scoats-wrf/{Archive}/netcdf_")
    wrf_links <- map(seq(year(dateRng[1]), year(dateRng[2]), by=1), 
                     ~glue("{thredds_head}{.x}{catalog}") |>
                       read_html() |> html_nodes("a") |> html_attr("href") |> 
                       grep("netcdf_20.*/wrf_", x=_, value=T) |>
                       str_split_fixed(glue("{Archive}/"), 2) |> 
                       magrittr::extract(,2)) |>
      do.call('c', args=_)
    wrf_i <- tibble(fname=wrf_links)
    wrf_base <- glue("{thredd_base}dodsC/scoats-wrf/{Archive}")
  } else {
    wrf_i <- tibble(fname=dir(wrf.dir, ".nc$", recursive=T))
    wrf_base <- wrf.dir
  }
  wrf_i <- wrf_i |>
    mutate(res=str_sub(fname, -6-forecast, -4-forecast),
           day_0=str_sub(fname, 23+forecast, 24+forecast),
           month_0=str_sub(fname, 21+forecast, 22+forecast),
           year_0=str_sub(fname, 17+forecast, 20+forecast),
           day_1=str_sub(fname, 28+forecast, 29+forecast),
           month_1=str_sub(fname, 26+forecast, 27+forecast),
           year_1=if_else(month_0=="12" & month_1=="01",
                          as.character(as.numeric(year_0) + 1),
                          year_0),
           date_0=ymd(paste0(year_0, month_0, day_0)),
           date_1=ymd(paste0(year_1, month_1, day_1))) |>
    filter(date_0 >= dateRng[1]-nDays_buffer, 
           date_1 <= dateRng[2]+nDays_buffer)
  wrf_dates <- unique(wrf_i$date_0)
  
  for(i in 1:length(wrf_dates)) {
    wrf_i.i <- filter(wrf_i, date_0==wrf_dates[i])
    nc_f.i <- map(wrf_i.i$fname, ~glue("{wrf_base}/{.x}")) |> 
      set_names(wrf_i.i$res)
    nc.ls <- map(nc_f.i, nc_open)
    time.ls <- map(nc.ls, 
                   ~tibble(Times=ncvar_get(.x, "Times")) |>
                     mutate(Time.dt=as_datetime(str_replace(Times, "_", " "))))
    var.ls <- map2_dfr(nc.ls, time.ls,
                       ~expand_grid(time=.y$Time.dt,
                                    lat_i=1:(.x$dim$south_north$len),
                                    lon_i=1:(.x$dim$west_east$len)) |>
                         mutate(date=date(time),
                                U=c(ncvar_get(.x, "U10")),
                                V=c(ncvar_get(.x, "V10")),
                                Shortwave=c(ncvar_get(.x, "Shortwave")),
                                Precipitation=c(ncvar_get(.x, "Precipitation")),
                                sst=c(ncvar_get(.x, "sst"))) |>
                         group_by(time) |>
                         mutate(i=row_number()) |>
                         ungroup(), 
                       .id="res") |> 
      group_by(res, date) |>
      group_split()
    for(j in seq_along(var.ls)) {
      j.fname <- glue("{out.dir}/wrf_{var.ls[[j]]$date[1]}_{var.ls[[j]]$res[1]}.rds")
      j.fnameF <- glue("{out.dir}/wrfF_{var.ls[[j]]$date[1]}_{var.ls[[j]]$res[1]}.rds")
      if(!file.exists(j.fname)) {
        var.ls[[j]] |> 
          mutate(sst=if_else(sst > -100, sst, NA_real_)) |>
          group_by(i, lat_i, lon_i, date) |>
          summarise(U_mn=mean(U),
                    V_mn=mean(V),
                    UV_mn=mean(sqrt(U^2 + V^2)),
                    Shortwave=sum(Shortwave),
                    Precip=sum(Precipitation),
                    sst=mean(sst, na.rm=T)) |>
          ungroup() |>
          rename(U=U_mn, V=V_mn, UV=UV_mn) |>
          saveRDS(ifelse(forecast, j.fnameF, j.fname))
      }
    }
    # save domain extents IF all three domains are represented
    if(length(nc_f.i)==3 & all(c("d01", "d02", "d03") %in% names(nc_f.i))) {
      nest_WRF_domains(nc.ls) |>
        walk(~saveRDS(.x, glue("{out.dir}/wrfDomains_{wrf_dates[i]}_{.x$res[1]}.rds")))
    }
    walk(nc.ls, nc_close)
  }
}






#' Identify preserved points within nested WRF domains 
#'
#' @param nc.ls 
#'
#' @return
#' @export
#'
#' @examples
nest_WRF_domains <- function(nc.ls) {
  library(tidyverse); library(ncdf4); library(sf)
  nDomain <- length(nc.ls)
  lon <- map(nc.ls, ~ncvar_get(.x, "XLONG"))
  lat <- map(nc.ls, ~ncvar_get(.x, "XLAT"))
  elev <- map(nc.ls, ~ncvar_get(.x, "HGT"))
  coord.rows <- map(lon, ~matrix(1:nrow(.x), nrow=nrow(.x), ncol=ncol(.x)))
  coord.cols <- map(lon, ~matrix(1:ncol(.x), nrow=nrow(.x), ncol=ncol(.x), byrow=T))
  coord.wrf <- map(1:nDomain, 
                   ~tibble(lon=c(lon[[.x]]),
                           lat=c(lat[[.x]]),
                           elev=c(elev[[.x]]),
                           row=c(coord.rows[[.x]]),
                           col=c(coord.cols[[.x]])) |>
                     mutate(i=row_number(),
                            res=names(nc.ls)[.x]) |>
                     st_as_sf(coords=c("lon", "lat"), crs=4326))
  hull.wrf <- map(1:nDomain, 
                  ~st_convex_hull(st_union(coord.wrf[[.x]])) |>
                    st_as_sf() |> mutate(res=names(nc.ls)[.x]))
  merge.wrf <- map(coord.wrf, ~.x |> mutate(in_d02=0, in_d03=0))
  if(nDomain==3) {
    merge.wrf <- map(merge.wrf, 
                     ~.x %>%
                       mutate(in_d02=as.numeric(st_within(., hull.wrf[[2]], sparse=F)[,1]),
                              in_d03=as.numeric(st_within(., hull.wrf[[3]], sparse=F)[,1])))
  }
  if(nDomain==2) {
    merge.wrf <- map(merge.wrf, 
                     ~.x %>%
                       mutate(in_d02=as.numeric(st_within(., hull.wrf[[2]], sparse=F)[,1])))
  }
  merge.wrf <- map(merge.wrf,
                   ~.x %>%
                     mutate(lon=st_coordinates(.)[,1],
                            lat=st_coordinates(.)[,2]) |>
                     st_drop_geometry() |>
                     filter(res=="d03" | (res=="d02" & !in_d03) | (res=="d01" & !in_d02)) |>
                     select(-in_d02, -in_d03))
  return(merge.wrf)
}





#' Filter WRF data to include only preserved points in each domain
#'
#' @param domain 
#' @param wrf.out 
#' @param v2_start 
#'
#' @return
#' @export
#'
#' @examples
subset_WRF <- function(domain, wrf.out, v2_start=NULL, refreshStart=NULL) {
  f.domain <- dirf(wrf.out, glue("wrfDomains_.*{domain}.rds"))
  f.wrf <- dirf(wrf.out, glue("wrfF?_.*{domain}.rds"))
  if(is.null(v2_start)) {
    domain.ls <- map(f.domain, readRDS)
    i_chg <- c(1, which(map_lgl(1:(length(domain.ls)-1), 
                                ~!(identical(domain.ls[[.x]]$lon, domain.ls[[.x+1]]$lon) &
                                     identical(domain.ls[[.x]]$lat, domain.ls[[.x+1]]$lat))))+1)
    domain.ls <- domain.ls[i_chg]
  } else {
    i_chg <- c(1, grep(v2_start, f.domain))
    domain.ls <- map(f.domain[i_chg], readRDS)
  }
  v_dateRng <- tibble(v=seq_along(i_chg),
                      start=str_split_fixed(str_split_fixed(f.domain[i_chg], "Domains_", 2)[,2], "_d0", 2)[,1]) |>
    mutate(start=ymd(start), 
           end=lead(start, default=ymd("3000-01-01")))
  v_dateRng$start[1] <- "2013-01-01"
  
  iwalk(domain.ls, 
        ~.x |> mutate(version=.y) |>
          select(-row, -col) |>
          saveRDS(glue("{wrf.out}/domain_{domain}_{.y}.rds")))
  
  wrf.ls <- vector("list", length(f.wrf))
  for(i in seq_along(f.wrf)) {
    date_i <- str_split_fixed(str_split_fixed(f.wrf[i], "/wrfF?_", 2)[,2], "_d0", 2)[,1]
    if(is.null(refreshStart) || date_i >= refreshStart) {
      v_i <- which(date_i >= v_dateRng$start & date_i < v_dateRng$end)
      wrf.ls[[i]] <- readRDS(f.wrf[i]) |>
        select(-lon_i, -lat_i) |>
        right_join(domain.ls[[v_i]] |> select(-row, -col), by="i") |>
        mutate(version=ifelse(is.null(v2_start), v_i, 1 + (date_i >= v2_start))) 
    }
  }
  wrf.ls <- do.call('rbind', wrf.ls)
  return(wrf.ls)
}





aggregate_WRF <- function(wrf.out, v2_start=ymd("2019-04-01"), refreshStart=NULL) {
  wrf.df <- subset_WRF("d01", wrf.out, v2_start=ymd("2019-04-01"), refreshStart) |>
    bind_rows(subset_WRF("d02", wrf.out, v2_start=ymd("2019-04-01"), refreshStart)) |>
    bind_rows(subset_WRF("d03", wrf.out, v2_start=ymd("2019-04-01"), refreshStart)) |>
    filter(!is.na(date)) |>
    arrange(date, res, i) |>
    group_by(date) |>
    mutate(wrf_id=row_number()) |>
    ungroup() |>
    mutate(across(where(is.numeric), ~if_else(.x > 1e30, NA, .x))) |>
    mutate(yday=yday(date)) |>
    group_by(wrf_id, version, yday) |>
    mutate(across(where(is.numeric), zoo::na.aggregate)) |>
    ungroup() |>
    mutate(Shortwave=log1p(Shortwave),
           Precip=log1p(pmax(Precip, 0)*3600*24*1000), # m/s to mm/day
           UV=log1p(UV))
  return(wrf.df)
}






#' Extract CMEMS or WRF data to site point locations
#'
#' @param site.df 
#' @param env_vars 
#' @param env.df 
#' @param id_env 
#' @param site.v 
#'
#' @return
#' @export
#'
#' @examples
extract_env_pts <- function(site.df, env_vars, env.df, id_env, site.v) {
  library(tidyverse); library(zoo)
  
  env.site <- env.df |>
    group_by(version) |>
    group_split() |>
    map2_dfr(.x=_, .y=site.v, ~.x |> filter({{id_env}} %in% site.df[[.y]])) |>
    arrange({{id_env}}, date) |>
    group_by({{id_env}}) |>
    mutate(across(any_of(env_vars), 
                  ~rollmeanr(.x, k=7, na.pad=T, align="right"),
                  .names="{.col}Wk")) |>
    mutate(across(any_of(paste0(env_vars, "Wk")),
                  ~.x - lag(.x),
                  .names="{.col}Delta")) |>
    # ungroup() |>
    # mutate(yday=yday(date)) |>
    # group_by({{id_env}}) |>
    # mutate(across(any_of(env_vars),
    #               ~detrend_loess(yday, .x, span=0.3), 
    #               .names="{.col}Dt")) |>
    ungroup()
  env.site <- env.site |> 
    select({{id_env}}, version, date, 
           any_of(paste0(env_vars, "Wk")),
           any_of(paste0(env_vars, "WkDelta")),
           any_of(paste0(env_vars, "Dt")))
  return(env.site)
}







#' Extract CMEMS or WRF data within site buffer quadrants
#'
#' @param site.buffer 
#' @param vars 
#' @param env.df 
#' @param id_env 
#'
#' @return
#' @export
#'
#' @examples
extract_env_buffers <- function(site.buffer, vars, env.df, id_env) {
  
  library(tidyverse); library(zoo)
  
  env.df <- env.df |> arrange(date, pick(any_of(id_env)))
  env.buffer <- expand_grid(siteid=unique(site.buffer$siteid),
                            quadrant=unique(site.buffer$quadrant),
                            date=unique(env.df$date)) |>
    mutate(v=1 + ((date >= "2019-04-01")*(length(id_env)>1))) |>
    bind_cols(as_tibble(setNames(map(vars$all, ~NA_real_), vars$all)))
  
  env_id.ls <- map(id_env, ~map(site.buffer[[.x]], ~.x))
  env.df$date_id <- match(env.df$date, unique(env.df$date)) 
  env_dates.ls <- split(1:nrow(env.df), env.df$date_id)
  
  ij <- 1
  startTime <- Sys.time()
  if(is.null(vars$sea)) {
    for(i in 1:nrow(site.buffer)) {
      for(j in 1:length(env_dates.ls)) {
        if(length(env_id.ls[[env.buffer$v[ij]]][[i]]) > 0) {
          for(k in vars$all) {
            env.buffer[ij,k] <- mean(env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][[k]])
          }
        }
        ij <- ij+1
      }
      if(i %% 50 == 0) {
        cat(i, "of", nrow(site.buffer), "--", 
            round(as.numeric(difftime(Sys.time(), startTime, units="min")), 2), 
            "minutes \n")
      }
    }
  } else {
    for(i in 1:nrow(site.buffer)) {
      for(j in 1:length(env_dates.ls)) {
        if(length(env_id.ls[[env.buffer$v[ij]]][[i]]) > 0) {
          for(k in vars$sea) {
            env.buffer[ij,k] <- mean(env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][[k]])
          }
          sst_ij <- env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][["sst"]]
          elev_ij <- env.df[env_dates.ls[[j]],][env_id.ls[[env.buffer$v[ij]]][[i]],][["elev"]]
          env.buffer[ij,"sst"] <- mean(sst_ij[elev_ij==0], na.rm=T) 
        }
        ij <- ij+1
      }
      if(i %% 50 == 0) {
        cat(i, "of", nrow(site.buffer), "--", 
            round(as.numeric(difftime(Sys.time(), startTime, units="min")), 2), 
            "minutes \n")
      }
    }
  }
  
  env.buffer <- env.buffer |> 
    group_by(siteid, quadrant) |>
    mutate(across(any_of(vars$all), 
                  ~rollmean(.x, k=7, na.pad=T, align="right"),
                  .names="{.col}AvgWk")) |>
    mutate(across(any_of(paste0(vars$all, "AvgWk")),
                  ~.x - lag(.x),
                  .names="{.col}Delta")) |>
    # ungroup() |>
    # mutate(yday=yday(date)) |>
    # group_by(siteid, quadrant) |>
    # mutate(across(any_of(vars$all),
    #               ~detrend_loess(yday, .x, span=0.3), 
    #               .names="{.col}AvgDt")) |>
    ungroup() |> 
    select(siteid, quadrant, date, 
           any_of(paste0(vars$all, "AvgWk")),
           any_of(paste0(vars$all, "AvgWkDelta")),
           any_of(paste0(vars$all, "AvgDt"))) |>
    group_by(siteid, date) |>
    mutate(across(where(is.numeric), na.aggregate)) |>
    ungroup()
  return(env.buffer)
}








#' Find id of environmental layer points nearest to each site point location
#'
#' @param site.df 
#' @param env.sf 
#' @param id_env 
#'
#' @return
#' @export
#'
#' @examples
find_nearest_feature_id <- function(site.df, env.sf, id_env) {
  site.df |>
    st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) |>
    st_transform(4326) %>%
    mutate(new_id=st_nearest_feature(., env.sf)) |>
    rename_with(~id_env, new_id) |>
    st_drop_geometry()
}




#' Find ids of environmental layer points within each site buffer
#'
#' @param site.sf 
#' @param env.sf 
#' @param id_env 
#'
#' @return
#' @export
#'
#' @examples
find_buffer_intersect_ids <- function(site.sf, env.sf, id_env) {
  site.sf |>
    select(siteid, quadrant, geom) |>
    st_transform(4326) |>
    st_make_valid() %>%
    mutate(new_id=st_intersects(., env.sf)) |>
    rename_with(~id_env, new_id) |>
    st_drop_geometry()
}





#' Find pairwise shortest in-ocean paths between sites 
#'
#' @param ocean.path 
#' @param site.df 
#' @param transMx.path 
#' @param recalc_transMx 
#'
#' @return
#' @export
#'
#' @examples
get_shortestPaths <- function(ocean.path, site.df, site_savePath=NULL) {
  library(raster); library(gdistance); 
  library(terra); library(spaths)
  library(tidyverse); library(sf); library(glue);
  
  # adapted from:
  # https://agrdatasci.github.io/gdistance/reference/index.html
  
  # load ocean raster and calculate transition matrix
  mesh.r <- raster(ocean.path)
  crs(mesh.r) <- CRS("+init=epsg:27700")
  
  # locate sites in mesh
  set_ll_warn(TRUE)
  site.spdf <- SpatialPointsDataFrame(site.df[,c("lon", "lat")],
                                      data=site.df[,"siteid"], 
                                      proj4string=CRS("+init=epsg:27700")) |>
    points2nearestcell(mesh.r) |>
    as.data.frame()
  site.df_new <- site.df |> dplyr::select(siteid, sin) |> left_join(site.spdf)
  if(!is.null(site_savePath)) {
    saveRDS(site.df_new, site_savePath)
  }
  
  # Pairwise each i to all others within site.df
  mesh.r <- rast(ocean.path)
  crs(mesh.r) <- crs("+init=epsg:27700")
  site.sf <- site.df_new |>
    st_as_sf(coords=c("lon", "lat"), crs=27700) |>
    st_crop(mesh.r)
  out_paths <- spaths::shortest_paths(mesh.r, 
                              site.sf, 
                              output="distances")
  # dist.df <- map_dfr(1:nrow(site.spdf),
  #                    ~shortestPath(mesh.tmx,
  #                                  as.matrix(site.spdf[.x, c("lon", "lat")]),
  #                                  as.matrix(site.spdf[-.x, c("lon", "lat")]),
  #                                  output="SpatialLines") |>
  #                      st_as_sf() |>
  #                      mutate(origin=site.spdf$siteid[.x],
  #                             destination=site.spdf$siteid[-.x],
  #                             distance=st_length(.)) |>
  #                      st_drop_geometry()) 
  
  return(list(site.df=site.df_new, dist.df=out_paths))
}





#' Identify whether point locations fall within valid raster cells
#'
#' @param locs 
#' @param ras 
#' @param layer 
#'
#' @return
#' @export
#'
#' @examples
point_in_cell <- function (locs, ras, layer=1) {
  # copied from rSDM since not available for newer R versions
  if (!isTRUE(raster::compareCRS(locs, ras))) {
    stop("Coordinate data and raster object must have the same projection. Check their CRS or proj4string")
  }
  else {
    if (nlayers(ras) > 1) 
      ras <- raster(ras, layer)
    rasvals <- raster::extract(ras, locs)
    missing <- is.na(rasvals)
    missing
  }
}



#' Shift out-of-bounds pointsto nearest cell
#'
#' @param locs 
#' @param ras 
#' @param layer 
#' @param move 
#' @param distance 
#' @param showchanges 
#'
#' @return
#' @export
#'
#' @examples
points2nearestcell <- function (locs=NULL, ras=NULL, layer=1, 
                                move=T, distance=NULL, showchanges=T) {
  # copied from rSDM since not available for newer R versions
  miss <- point_in_cell(locs, ras, layer)
  if (sum(miss) > 0) {
    coord.miss <- sp::coordinates(locs[miss, ])
    if (nlayers(ras) > 1) 
      ras <- raster::raster(ras, layer)
    cells.notNA <- raster::rasterToPoints(ras, spatial = TRUE)
    coord.ras <- sp::coordinates(cells.notNA)
    cell.id <- factor(seq_len(nrow(coord.ras)))
    nearest.cell <- class::knn1(coord.ras, coord.miss, cell.id)
    new.coords <- matrix(coord.ras[nearest.cell, ], ncol = 2)
    colnames(new.coords) <- c("longitude_new", "latitude_new")
    if (!is.null(distance)) {
      distances <- raster::pointDistance(coord.miss, new.coords, 
                                         lonlat = raster::isLonLat(locs))
      x <- ifelse(distances < distance, new.coords[, 1], 
                  coord.miss[, 1])
      y <- ifelse(distances < distance, new.coords[, 2], 
                  coord.miss[, 2])
      new.coords <- cbind(longitude_new = x, latitude_new = y)
    }
    if (isTRUE(move)) {
      locs@coords[miss, ] <- new.coords
    }
    if (isTRUE(showchanges)) {
      coords <- data.frame(coord.miss, new.coords)
      distances <- round(raster::pointDistance(coord.miss, 
                                               new.coords, lonlat = raster::isLonLat(locs)))
      moved <- apply(coords, 1, function(x) {
        !isTRUE(identical(x[1], x[3]) & identical(x[2], 
                                                  x[4]))
      })
      coords <- cbind(coords, distances, moved)
      print(coords)
      message(sum(moved), " out of ", nrow(coords), 
              " points have been moved.")
    }
  }
  else message("All points fall within a raster cell")
  return(locs)
}





#' Extract fetch for each site point location
#'
#' @param site.df 
#' @param fetch.path 
#'
#' @return
#' @export
#'
#' @examples
get_fetch <- function(site.df, fetch.path) {
  library(tidyverse); library(sf)
  site.df |>
    st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) %>%
    mutate(fetch=raster::extract(raster(fetch.path), ., 
                                 small=T, buffer=1e2, fun=mean)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), ., 
                                         small=T, buffer=5e2, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), ., 
                                         small=T, buffer=1e3, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), ., 
                                         small=T, buffer=5e3, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), ., 
                                         small=T, buffer=10e3, fun=mean),
                         fetch)) %>%
    mutate(fetch=if_else(is.na(fetch),
                         raster::extract(raster(fetch.path), ., 
                                         small=T, buffer=50e3, fun=mean),
                         fetch)) %>%
    st_drop_geometry()
}





#' Identify most-open bearing for each site point location
#'
#' @param site.df 
#' @param coast.path 
#' @param buffer 
#' @param nDir 
#'
#' @return
#' @export
#'
#' @examples
get_openBearing <- function(site.df, coast.path, buffer=10e3, nDir=120) {
  library(tidyverse); library(sf)
  
  hub.df <- site.df |> 
    select(siteid, lon, lat) |>
    st_as_sf(coords=c("lon", "lat"), crs=27700)
  spoke.df <- hub.df |>
    st_buffer(dist=buffer, nQuadSegs=nDir/4) |>
    st_cast("POINT") |>
    group_by(siteid) |>
    mutate(spoke.id=row_number()) |>
    filter(spoke.id <= nDir) 
  hubRep.df <- full_join(hub.df, spoke.df |> st_drop_geometry())
  coords <- cbind(st_coordinates(hubRep.df), st_coordinates(spoke.df))
  spoke.lines <- lapply(1:nrow(coords),
                        function(i){
                          st_linestring(matrix(coords[i,], ncol=2, byrow=TRUE))
                        }) |>
    st_sfc() |> st_as_sf() |> st_set_crs(27700) |>
    rename(geometry=x) |>
    mutate(siteid=spoke.df$siteid,
           spoke.id=spoke.df$spoke.id)
  spoke.mesh <- st_intersection(spoke.lines, st_read(coast.path)) |>
    st_cast("LINESTRING") %>%
    mutate(len=round(as.numeric(st_length(.)))) |>
    group_by(siteid) |>
    filter(len==max(len)) |>
    ungroup() |>
    st_cast("POINT")
  bearings <- as.numeric(lwgeom::st_geod_azimuth(st_transform(spoke.mesh, 4326)))
  bearing.df <- spoke.mesh[1:nrow(spoke.mesh) %% 2 == 0,] |>
    st_drop_geometry() |>
    mutate(bearing=bearings[1:length(bearings) %% 2 == 1]) |>
    group_by(siteid) |>
    summarise(openBearing=median(bearing)) |>
    ungroup()
  return(full_join(site.df, bearing.df, by="siteid"))
}





#' Split circular buffer into NSEW quadrants
#'
#' @param sf 
#' @param radius 
#'
#' @return
#' @export
#'
#' @examples
split_to_NSEW <- function(sf, radius=110e3) {
  library(tidyverse); library(sf); library(lwgeom)
  hub.df <- sf |> 
    st_centroid() |>
    select(siteid, geometry)
  spoke.df <- hub.df |>
    st_buffer(dist=radius, nQuadSegs=2) |>
    st_cast("POINT") |>
    group_by(siteid) |>
    mutate(spoke.id=row_number()) |>
    filter(spoke.id %% 2 == 0) |>
    mutate(side=c("start", "start", "end", "end")) |>
    ungroup()
  coords <- cbind(st_coordinates(filter(spoke.df, side=="start")),
                  st_coordinates(filter(spoke.df, side=="end")))
  spoke.lines <- map(1:nrow(coords),
                     ~st_linestring(matrix(coords[.x,],ncol=2, byrow=T))) |>
    st_sfc() |> st_as_sf() |> st_set_crs(27700) |>
    rename(geometry=x) |>
    mutate(siteid=filter(spoke.df, side=="start")$siteid)
  sf.quad <- sf$siteid |>
    map_dfr(~st_split(filter(sf, siteid==.x), filter(spoke.lines, siteid==.x)) |>
              st_collection_extract() |>
              mutate(quadrant=c("E", "S", "W", "N"))) 
  
  return(sf.quad)
}





#' Lag multiple variables at once
#'
#' From https://stackoverflow.com/questions/55814028/multiple-lags-with-dplyr
#'
#' @param data Dataframe
#' @param ... Unquoted variable names to lag
#' @param n Number of lags
#'
#' @return
#' @export
#'
#' @examples
get_lags <- function(data, ..., n=2){
  library(tidyverse); library(rlang)
  variable <- enquos(...)
  
  indices <- seq_len(n)
  combos <- tidyr::crossing(indices, var=as.list(variable))
  
  quosures <- map2(combos$indices, combos$var,
                   ~quo(lag(!!.y, !!.x)) ) |>
    set_names(paste0(map_chr(combos$var, quo_text), combos$indices))
  mutate(data, !!!quosures )
}





#' Identify traffic light and alert status based on density/concentration
#'
#' @param y.df 
#' @param N 
#' @param tl_i 
#'
#' @return
#' @export
#'
#' @examples
get_trafficLights <- function(y.df, N, tl_i) {
  library(tidyverse)
  y.df |>
    rowwise() |>
    mutate(tl=tl_i$tl[max(which(y==tl_i$abbr & {{N}} >= tl_i$min_ge))],
           alert=tl_i$A[max(which(y==tl_i$abbr & {{N}} >= tl_i$min_ge))]) |>
    ungroup()
}





#' Calculate local and regional autoregressive terms for HAB and toxin observations
#'
#' @param yRaw.df 
#' @param y_i 
#' @param tl_i 
#' @param site.100km 
#'
#' @return
#' @export
#'
#' @examples
calc_y_features <- function(yRaw.df, y_i, tl_i, dist.df=NULL, forecastStart="3000-01-01") {
  y.ls <- yRaw.df |> 
    group_by(siteid, date) |>
    slice_head(n=1) |>
    ungroup() |>
    pivot_longer(any_of(y_i$abbr), names_to="y", values_to="N") |>
    filter((!is.na(N)) | (date >= forecastStart)) |>
    mutate(lnN=log1p(N)) |>
    get_trafficLights(N, tl_i) |>
    arrange(y, siteid, date) |>
    group_by(y, siteid) |>
    get_lags(lnN, alert, date, n=2) |>
    ungroup() |>
    mutate(year=year(date),
           lnNWt1=lnN1/log1p(as.numeric(date-date1)), 
           lnNWt2=lnN2/log1p(as.numeric(date-date2)),
           lnNAvg1=0, lnNAvg2=0, prAlertAvg1=0, prAlertAvg2=0, 
           lnNAvgPrevYr=0, prAlertAvgPrevYr=0) |>
    group_by(y) |>
    group_split()
  y_new <- vector("list", length(y.ls))
  
  for(i in seq_along(y.ls)) {
    y.df_i <- y.ls[[i]] |> select(siteid, date, lnN1, lnN2, alert1, alert2)
    yYr_i <- y.ls[[i]] |> 
      group_by(siteid, year) |>
      summarise(lnN_yr=mean(lnN, na.rm=T),
                prAlert_yr=mean(alert=="A1", na.rm=T)) |>
      ungroup() |>
      mutate(year_matchObs=year+1)
    y_new[[i]] <- y.ls[[i]] |>
      left_join(yYr_i |> select(year_matchObs, siteid, lnN_yr, prAlert_yr), 
                by=c("year"="year_matchObs", "siteid"="siteid")) |>
      rename(lnNPrevYr=lnN_yr, prAlertPrevYr=prAlert_yr)
    
    for(j in 1:nrow(y.ls[[i]])) {
      site_j <- y.df_i$siteid[j]
      date_j <- y.df_i$date[j]
      yr_j <- year(date_j)
      wk.df <- y.df_i |>
        filter(siteid %in% dist.df$dest_c[dist.df$origins==site_j][[1]] &
                 date <= date_j & 
                 date > date_j-7) 
      yr.df <- yYr_i |>
        filter(siteid %in% dist.df$dest_c[dist.df$origins==site_j][[1]] &
                 year == yr_j - 1)
      yr.site_j <- yYr_i |> filter(siteid==site_j & year==yr_j-1)
      y_new[[i]]$lnNAvg1[j] <- mean(wk.df$lnN1, na.rm=T)
      y_new[[i]]$lnNAvg2[j] <- mean(wk.df$lnN2, na.rm=T)
      y_new[[i]]$prAlertAvg1[j] <- mean(wk.df$alert1 == "A1", na.rm=T)
      y_new[[i]]$prAlertAvg2[j] <- mean(wk.df$alert2 == "A1", na.rm=T)
      y_new[[i]]$lnNAvgPrevYr[j] <- mean(yr.df$lnN_yr, na.rm=T)
      y_new[[i]]$prAlertAvgPrevYr[j] <- mean(yr.df$prAlert_yr, na.rm=T)
      if(j %% 1000 == 0) {cat(i, ":", j, "of", nrow(y_new[[i]]), "\n")}
    }
  }
  y.df <- do.call('rbind', y_new)
  return(y.df)
}




#' Summarise HAB states for each toxin observation
#'
#' @param site_tox.sf 
#' @param site_hab.sf 
#' @param tox.obs 
#' @param hab.df 
#'
#' @return
#' @export
#'
#' @examples
summarise_hab_states <- function(site_tox.sf, site_hab.sf, tox.obs, hab.df) {
  library(tidyverse); library(sf)
  
  # identify hab sites within each toxin site buffer
  hab_ids <- site_tox.sf |>
    select(siteid, geom) %>%
    mutate(hab_id=st_intersects(., site_hab.sf)) |>
    st_drop_geometry()
  
  # match hab observations to toxin observations
  hab.df <- hab.df |>
    mutate(prA=as.numeric(alert=="A1")) |>
    select(siteid, date, y, lnN, prA) |> 
    rename(hab_id=siteid, lnNAvg=lnN) |>
    pivot_wider(names_from=y, values_from=c(lnNAvg, prA), names_glue="{y}{.value}")
  hab_y_names <- names(hab.df)[-(1:2)]
  habSums <- vector("list", nrow(tox.obs))
  for(i in 1:nrow(tox.obs)) {
    hab_sites <- filter(hab_ids, siteid==tox.obs$siteid[i])$hab_id[[1]]
    habSums[[i]] <- hab.df |> 
      filter(hab_id %in% hab_sites & 
               date <= tox.obs$date[i] - 7*1 &
               date >= tox.obs$date[i] - 7*5) |>
      summarise(across(all_of(hab_y_names), ~mean(.x, na.rm=T)))
  }
  out.df <- tox.obs |> bind_cols(do.call('rbind', habSums))
  return(out.df)
}





#' Detrend observations using a loess smoother
#'
#' @param x 
#' @param y 
#' @param span 
#' @param robust 
#'
#' @return
#' @export
#'
#' @examples
detrend_loess <- function (x, y, span=0.75, robust=TRUE) {
  # modified from astsa::trend
  if(sum(!is.na(y)) < 10) {
    return(y)
  }
  if(length(y) < 10) {
    return(y)
  }
  fam = ifelse(robust, "symmetric", "gaussian")
  lo = stats::predict(stats::loess(y ~ x, span=span, family=fam), se = F)
  return(c(y - lo))
}






#' Load datasets for compilation
#'
#' @param sub.dir 
#' @param target 
#'
#' @return
#' @export
#'
#' @examples
load_datasets <- function(sub.dir, target) {
  d.ls <- list(
    site=readRDS(glue("data/site_{target}_df.rds")),
    obs=readRDS(glue("data/{sub.dir}/{target}_obs.rds")),
    cmems.pt=readRDS(glue("data/{sub.dir}/cmems_sitePt_{target}.rds")),
    cmems.buf=readRDS(glue("data/{sub.dir}/cmems_siteBufferNSEW_{target}.rds")),
    wrf.pt=readRDS(glue("data/{sub.dir}/wrf_sitePt_{target}.rds")),
    wrf.buf=readRDS(glue("data/{sub.dir}/wrf_siteBufferNSEW_{target}.rds")),
    fsa=readRDS(glue("data/{sub.dir}/fsa_df.rds")),
    cefas=readRDS(glue("data/{sub.dir}/cefas_df.rds"))
  )
  if(target=="tox") {
    d.ls$habAvg <- readRDS(glue("data/{sub.dir}/tox_habAvg.rds"))
  }
  yday_env <- list(
    cmems.pt=readRDS(glue("data/1_current/ydayAvg_cmems_sitePt_{target}.rds")) |> 
      filter(version==max(version)) |> select(-version),
    cmems.buf=readRDS(glue("data/1_current/ydayAvg_cmems_siteBufferNSEW_{target}.rds")) |>
      pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
    wrf.pt=readRDS(glue("data/1_current/ydayAvg_wrf_sitePt_{target}.rds")) |> 
      filter(version==max(version)) |> select(-version),
    wrf.buf=readRDS(glue("data/1_current/ydayAvg_wrf_siteBufferNSEW_{target}.rds")) |>
      pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir")
  )
  d.ls$compiled <- d.ls$site |> select(-sin) |>
    right_join(d.ls$obs, by="siteid", multiple="all") |>
    left_join(d.ls$cmems.pt |> select(-ends_with("Dt")), by=c("cmems_id", "date")) |>
    left_join(d.ls$cmems.buf |> select(-ends_with("Dt")) |>
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
              by=c("siteid", "date")) |>
    mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) |>
    select(-wrf_id.1, -wrf_id.2, -version) |>
    left_join(d.ls$wrf.pt |> select(-version, -ends_with("Dt")), by=c("wrf_id", "date")) |>
    left_join(d.ls$wrf.buf |> select(-ends_with("Dt")) |>
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
              by=c("siteid", "date")) |>
    mutate(year=year(date),
           yday=yday(date))
  for(i in seq_along(yday_env)) {
    env_id_col <- switch(i,
                         "1"=c("cmems_id", "yday"),
                         "2"=c("siteid", "yday"),
                         "3"=c("wrf_id", "yday"),
                         "4"=c("siteid", "yday"))
    varNames_i <- names(yday_env[[i]] |> select(-all_of(env_id_col)))
    for(j in varNames_i) {
      na_rows <- which(is.na(d.ls$compiled[[j]]))
      if(length(na_rows > 0)) {
        d_meta <- d.ls$compiled[na_rows, c("siteid", "cmems_id", "wrf_id", "yday")]
        d.ls$compiled[[j]][na_rows] <- left_join(d_meta |> 
                                                   select(all_of(env_id_col)),
                                                 yday_env[[i]] |> 
                                                   select(all_of(c(env_id_col, j))))[[j]]
      }
    }
  }
  if(target=="tox") {
    d.ls$compiled <- d.ls$compiled |>
      left_join(d.ls$habAvg |> select(-date, -siteid), by="obsid")
  }
  
  return(d.ls)
}



