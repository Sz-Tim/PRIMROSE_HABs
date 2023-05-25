# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Functions







# Data preparation --------------------------------------------------------

read_and_clean_sites <- function(url_sites, dateStart) {
  url(url_sites) %>%
    readLines(warn=F) %>%
    fromJSON() %>% as_tibble %>%
    filter(east < 7e5,
           north < 125e4,
           !(east==0 & north==0),
           sin != "-99",
           !is.na(fromdate) & !is.na(todate),
           fromdate != todate,
           todate >= dateStart) %>%
    mutate(fromdate=date(fromdate), todate=date(todate)) %>%
    rowwise() %>%
    mutate(date=list(seq(fromdate, todate, by=1))) %>%
    ungroup %>%
    arrange(sin, fromdate) %>%
    select(sin, east, north, date) %>%
    unnest(date) %>%
    arrange(sin, date) %>%
    group_by(sin, date) %>%
    slice_head(n=1) %>%
    ungroup
}



read_and_clean_fsa <- function(url_fsa, hab_i, sites, dateStart="2016-01-01") {
  url(url_fsa) %>% 
    readLines(warn=F) %>%
    fromJSON() %>% as_tibble %>% 
    filter(!is.na(date_collected)) %>%
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) %>%
    filter(date >= dateStart) %>%
    mutate(across(any_of(hab_i$full), ~na_if(.x, -99))) %>%
    group_by(sin) %>% mutate(N=n()) %>% ungroup %>% filter(N > 2) %>%
    select(oid, sin, date, easting, northing, all_of(hab_i$full)) %>%
    left_join(sites, by=c("sin", "date")) %>%
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) %>%
    rename(obsid=oid) %>%
    group_by(sin) %>% mutate(lon=median(east), lat=median(north)) %>% ungroup %>%
    rename(all_of(setNames(hab_i$full, hab_i$abbr))) %>%
    select(obsid, lon, lat, sin, date, all_of(hab_i$abbr)) %>%
    arrange(sin, date) 
}


read_and_clean_cefas <- function(url_cefas, tox_i, sites, dateStart="2016-01-01") {
  url(url_cefas) %>% 
    readLines(warn=F) %>%
    fromJSON() %>% as_tibble %>% 
    filter(!is.na(date_collected) & sin != "-99") %>%
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) %>%
    filter(date >= dateStart) %>%
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
    rename(all_of(setNames(tox_i$full, tox_i$abbr))) %>%
    select(obsid, lon, lat, sin, date, all_of(tox_i$abbr)) %>%
    arrange(sin, date)
}


read_and_clean_fish <- function(url_mowi, url_ssf, fish_i, sites, dateStart="2016-01-01") {
  bind_rows(url(url_mowi) %>% 
              readLines(warn=F) %>%
              fromJSON() %>% as_tibble,
            url(url_ssf) %>% 
              readLines(warn=F) %>%
              fromJSON() %>% as_tibble) %>% 
    filter(!is.na(date_collected)) %>%
    mutate(datetime_collected=as_datetime(date_collected),
           date=date(datetime_collected)) %>%
    filter(date >= dateStart) %>%
    mutate(across(any_of(fish_i$full), ~na_if(.x, -99))) %>%
    group_by(sin) %>% mutate(N=n()) %>% ungroup %>% filter(N > 2) %>%
    select(oid, sin, date, easting, northing, all_of(fish_i$full)) %>%
    left_join(sites, by=c("sin", "date")) %>%
    mutate(east=if_else(is.na(east), easting, east),
           north=if_else(is.na(north), northing, north)) %>%
    rename(obsid=oid) %>%
    group_by(sin) %>% mutate(lon=median(east), lat=median(north)) %>% ungroup %>%
    rename(all_of(setNames(fish_i$full, fish_i$abbr))) %>%
    select(obsid, lon, lat, sin, date, all_of(fish_i$abbr)) %>%
    arrange(sin, date) %>%
    filter(sin!=0) %>%
    group_by(date, sin) %>% 
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
get_CMEMS <- function(userid, pw, i.df, bbox, nDays_buffer, dateRng, out.dir) {
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
                lon=c(nc.lon[lon_i])) %>%
      mutate(source=i.df$source[i],
             vals=c(nc.var)) %>%
      rename_with(~gsub("vals", i.df$var[i], .x)) %>%
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
                     ~glue("{thredds_head}{.x}{catalog}") %>%
                       read_html() %>% html_nodes("a") %>% html_attr("href") %>% 
                       grep("netcdf_20.*/wrf_", ., value=T) %>%
                       str_split_fixed(., glue("{Archive}/"), 2) %>% 
                       magrittr::extract(,2)) %>%
      do.call('c', .)
    wrf_i <- tibble(fname=wrf_links)
    wrf_base <- glue("{thredd_base}dodsC/scoats-wrf/{Archive}")
  } else {
    wrf_i <- tibble(fname=dir(wrf.dir, ".nc$", recursive=T))
    wrf_base <- wrf.dir
  }
  wrf_i <- wrf_i %>%
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
           date_1=ymd(paste0(year_1, month_1, day_1))) %>%
    filter(date_0 >= dateRng[1]-nDays_buffer, 
           date_1 <= dateRng[2]+nDays_buffer)
  wrf_dates <- unique(wrf_i$date_0)
  
  for(i in 1:length(wrf_dates)) {
    wrf_i.i <- filter(wrf_i, date_0==wrf_dates[i])
    nc_f.i <- map(wrf_i.i$fname, ~glue("{wrf_base}/{.x}")) %>% 
      set_names(wrf_i.i$res)
    nc.ls <- map(nc_f.i, nc_open)
    time.ls <- map(nc.ls, 
                   ~tibble(Times=ncvar_get(.x, "Times")) %>%
                     mutate(Time.dt=as_datetime(str_replace(Times, "_", " "))))
    var.ls <- map2_dfr(nc.ls, time.ls,
                       ~expand_grid(time=.y$Time.dt,
                                    lat_i=1:(.x$dim$south_north$len),
                                    lon_i=1:(.x$dim$west_east$len)) %>%
                         mutate(date=date(time),
                                U=c(ncvar_get(.x, "U10")),
                                V=c(ncvar_get(.x, "V10")),
                                Shortwave=c(ncvar_get(.x, "Shortwave")),
                                Precipitation=c(ncvar_get(.x, "Precipitation")),
                                sst=c(ncvar_get(.x, "sst"))) %>%
                         group_by(time) %>%
                         mutate(i=row_number()) %>%
                         ungroup, 
                       .id="res") %>% 
      group_by(res, date) %>%
      group_split()
    for(j in seq_along(var.ls)) {
      j.fname <- glue("{out.dir}/wrf_{var.ls[[j]]$date[1]}_{var.ls[[j]]$res[1]}.rds")
      j.fnameF <- glue("{out.dir}/wrfF_{var.ls[[j]]$date[1]}_{var.ls[[j]]$res[1]}.rds")
      if(!file.exists(j.fname)) {
        var.ls[[j]] %>% 
          mutate(sst=if_else(sst > -100, sst, NA_real_)) %>%
          group_by(i, lat_i, lon_i, date) %>%
          summarise(U_mn=mean(U),
                    V_mn=mean(V),
                    UV_mn=mean(sqrt(U^2 + V^2)),
                    Shortwave=sum(Shortwave),
                    Precip=sum(Precipitation),
                    sst=mean(sst)) %>%
          ungroup %>%
          rename(U=U_mn, V=V_mn, UV=UV_mn) %>%
          saveRDS(ifelse(forecast, j.fnameF, j.fname))
      }
    }
    # save domain extents IF all three domains are represented
    if(length(nc_f.i)==3 & all(c("d01", "d02", "d03") %in% names(nc_f.i))) {
      nest_WRF_domains(nc.ls) %>%
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
                           col=c(coord.cols[[.x]])) %>%
                     mutate(i=row_number(),
                            res=names(nc.ls)[.x]) %>%
                     st_as_sf(coords=c("lon", "lat"), crs=4326))
  hull.wrf <- map(1:nDomain, 
                   ~st_convex_hull(st_union(coord.wrf[[.x]])) %>%
                     st_as_sf %>% mutate(res=names(nc.ls)[.x]))
  merge.wrf <- map(coord.wrf, ~.x %>% mutate(in_d02=0, in_d03=0))
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
                            lat=st_coordinates(.)[,2]) %>%
                     st_drop_geometry() %>%
                     filter(res=="d03" | (res=="d02" & !in_d03) | (res=="d01" & !in_d02)) %>%
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
subset_WRF <- function(domain, wrf.out, v2_start=NULL) {
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
                      start=str_split_fixed(str_split_fixed(f.domain[i_chg], "Domains_", 2)[,2], "_d0", 2)[,1]) %>%
    mutate(start=ymd(start), 
           end=lead(start, default=ymd("3000-01-01")))
  v_dateRng$start[1] <- "2013-01-01"
  
  iwalk(domain.ls, 
        ~.x %>% mutate(version=.y) %>%
          select(-row, -col) %>%
          saveRDS(glue("{wrf.out}/domain_{domain}_{.y}.rds")))
  
  wrf.ls <- vector("list", length(f.wrf))
  for(i in seq_along(f.wrf)) {
    date_i <- str_split_fixed(str_split_fixed(f.wrf[i], "/wrfF?_", 2)[,2], "_d0", 2)[,1]
    v_i <- which(date_i >= v_dateRng$start & date_i < v_dateRng$end)
    wrf.ls[[i]] <- readRDS(f.wrf[i]) %>%
      select(-lon_i, -lat_i) %>%
      right_join(domain.ls[[v_i]] %>% select(-row, -col), by="i") %>%
      mutate(version=ifelse(is.null(v2_start), v_i, 1 + (date_i >= v2_start)))
  }
  wrf.ls <- do.call('rbind', wrf.ls)
  return(wrf.ls)
}





aggregate_WRF <- function(wrf.out, v2_start=ymd("2019-04-01")) {
  wrf.df <- subset_WRF("d01", wrf.out, v2_start=ymd("2019-04-01")) %>%
    bind_rows(subset_WRF("d02", wrf.out, v2_start=ymd("2019-04-01"))) %>%
    bind_rows(subset_WRF("d03", wrf.out, v2_start=ymd("2019-04-01"))) %>%
    filter(!is.na(date)) %>%
    arrange(date, res, i) %>%
    group_by(date) %>%
    mutate(wrf_id=row_number()) %>%
    ungroup %>%
    mutate(across(where(is.numeric), ~if_else(.x > 1e30, NA, .x))) %>%
    mutate(yday=yday(date)) %>%
    group_by(wrf_id, version, yday) %>%
    mutate(across(where(is.numeric), zoo::na.aggregate)) %>%
    ungroup %>%
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
  
  env.site <- env.df %>%
    group_by(version) %>%
    group_split() %>%
    map2_dfr(., site.v, ~.x %>% filter({{id_env}} %in% site.df[[.y]])) %>%
    arrange({{id_env}}, date) %>%
    group_by({{id_env}}) %>%
    mutate(across(any_of(env_vars), 
                  ~rollmeanr(.x, k=7, na.pad=T, align="right"),
                  .names="{.col}Wk")) %>%
    mutate(across(any_of(paste0(env_vars, "Wk")),
                  ~.x - lag(.x),
                  .names="{.col}Delta")) %>%
    ungroup %>%
    mutate(yday=yday(date)) %>%
    group_by({{id_env}}) %>%
    mutate(across(any_of(env_vars),
                  ~detrend_loess(yday, .x, span=0.3), 
                  .names="{.col}Dt")) %>%
    ungroup
  env.site <- env.site %>% 
    select({{id_env}}, version, date, 
           all_of(paste0(env_vars, "Wk")),
           all_of(paste0(env_vars, "WkDelta")),
           all_of(paste0(env_vars, "Dt")))
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
  
  env.df <- env.df %>% arrange(date, pick(any_of(id_env)))
  env.buffer <- expand_grid(siteid=unique(site.buffer$siteid),
                            quadrant=unique(site.buffer$quadrant),
                            date=unique(env.df$date)) %>%
    mutate(v=1 + ((date >= "2019-04-01")*(length(id_env)>1))) %>%
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
      cat(i, "of", nrow(site.buffer), "--", Sys.time()-startTime, "\n")
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
      cat(i, "of", nrow(site.buffer), "--", Sys.time()-startTime, "\n")
    }
  }
  
  env.buffer <- env.buffer %>% 
    group_by(siteid, quadrant) %>%
    mutate(across(any_of(vars$all), 
                  ~rollmean(.x, k=7, na.pad=T, align="right"),
                  .names="{.col}AvgWk")) %>%
    mutate(across(any_of(paste0(vars$all, "AvgWk")),
                  ~.x - lag(.x),
                  .names="{.col}Delta")) %>%
    ungroup %>%
    mutate(yday=yday(date)) %>%
    group_by(siteid, quadrant) %>%
    mutate(across(any_of(vars$all),
                  ~detrend_loess(yday, .x, span=0.3), 
                  .names="{.col}AvgDt")) %>%
    ungroup %>% 
    select(siteid, quadrant, date, 
           all_of(paste0(vars$all, "AvgWk")),
           all_of(paste0(vars$all, "AvgWkDelta")),
           all_of(paste0(vars$all, "AvgDt"))) %>%
    group_by(siteid, date) %>%
    mutate(across(where(is.numeric), na.aggregate)) %>%
    ungroup
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
  site.df %>%
    st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) %>%
    st_transform(4326) %>%
    mutate(new_id=st_nearest_feature(., env.sf)) %>%
    rename_with(~id_env, new_id) %>%
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
  site.sf %>%
    select(siteid, quadrant, geom) %>%
    st_transform(4326) %>%
    st_make_valid() %>%
    mutate(new_id=st_intersects(., env.sf)) %>%
    rename_with(~id_env, new_id) %>%
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
get_shortestPaths <- function(ocean.path, site.df, transMx.path=NULL, recalc_transMx=T, site_savePath=NULL) {
  library(tidyverse); library(sf); library(glue);
  library(raster); library(gdistance)
  
  # adapted from:
  # https://agrdatasci.github.io/gdistance/reference/index.html
  
  # load ocean raster and calculate transition matrix
  mesh.r <- raster(ocean.path)
  crs(mesh.r) <- CRS("+init=epsg:27700")
  if(is.null(transMx.path) | recalc_transMx) {
    mesh.tmx <- transition(mesh.r, mean, 16)
    mesh.tmx <- geoCorrection(mesh.tmx)
    saveRDS(mesh.tmx, "data/mesh_tmx.rds")
  } else {
    mesh.tmx <- readRDS(transMx.path)
  }
  
  # locate sites in mesh
  set_ll_warn(TRUE)
  site.spdf <- SpatialPointsDataFrame(site.df[,c("lon", "lat")],
                                      data=site.df[,"siteid"], 
                                      proj4string=CRS("+init=epsg:27700")) %>%
    points2nearestcell(., mesh.r) %>%
    as.data.frame
  site.df_new <- site.df %>% select(siteid, sin) %>% left_join(site.spdf)
  if(!is.null(site_savePath)) {
    saveRDS(site.df_new, site_savePath)
  }
  
  # Pairwise each i to all others within site.df
  dist.df <- map_dfr(1:nrow(site.spdf),
                     ~shortestPath(mesh.tmx,
                                   as.matrix(site.spdf[.x, c("lon", "lat")]),
                                   as.matrix(site.spdf[-.x, c("lon", "lat")]),
                                   output="SpatialLines") %>%
                       st_as_sf %>%
                       mutate(origin=site.spdf$siteid[.x],
                              destination=site.spdf$siteid[-.x],
                              distance=st_length(.)) %>%
                       st_drop_geometry) 
  
  return(list(site.df=site.df_new, dist.df=dist.df))
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
  site.df %>%
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
  
  hub.df <- site.df %>% 
    select(siteid, lon, lat) %>%
    st_as_sf(coords=c("lon", "lat"), crs=27700)
  spoke.df <- hub.df %>%
    st_buffer(dist=buffer, nQuadSegs=nDir/4) %>%
    st_cast("POINT") %>%
    group_by(siteid) %>%
    mutate(spoke.id=row_number()) %>%
    filter(spoke.id <= nDir) 
  hubRep.df <- full_join(hub.df, spoke.df %>% st_drop_geometry())
  coords <- cbind(st_coordinates(hubRep.df), st_coordinates(spoke.df))
  spoke.lines <- lapply(1:nrow(coords),
                        function(i){
                          st_linestring(matrix(coords[i,], ncol=2, byrow=TRUE))
                        }) %>%
    st_sfc() %>% st_as_sf() %>% st_set_crs(27700) %>%
    rename(geometry=x) %>%
    mutate(siteid=spoke.df$siteid,
           spoke.id=spoke.df$spoke.id)
  spoke.mesh <- st_intersection(spoke.lines, st_read(coast.path)) %>%
    st_cast("LINESTRING") %>%
    mutate(len=round(as.numeric(st_length(.)))) %>%
    group_by(siteid) %>%
    filter(len==max(len)) %>%
    ungroup %>%
    st_cast("POINT")
  bearings <- as.numeric(lwgeom::st_geod_azimuth(st_transform(spoke.mesh, 4326)))
  bearing.df <- spoke.mesh[1:nrow(spoke.mesh) %% 2 == 0,] %>%
    st_drop_geometry() %>%
    mutate(bearing=bearings[1:length(bearings) %% 2 == 1]) %>%
    group_by(siteid) %>%
    summarise(openBearing=median(bearing)) %>%
    ungroup
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
  hub.df <- sf %>% 
    st_centroid() %>%
    select(siteid, geometry)
  spoke.df <- hub.df %>%
    st_buffer(dist=radius, nQuadSegs=2) %>%
    st_cast("POINT") %>%
    group_by(siteid) %>%
    mutate(spoke.id=row_number()) %>%
    filter(spoke.id %% 2 == 0) %>%
    mutate(side=c("start", "start", "end", "end")) %>%
    ungroup
  coords <- cbind(st_coordinates(filter(spoke.df, side=="start")),
                  st_coordinates(filter(spoke.df, side=="end")))
  spoke.lines <- map(1:nrow(coords),
                     ~st_linestring(matrix(coords[.x,],ncol=2, byrow=T))) %>%
    st_sfc() %>% st_as_sf() %>% st_set_crs(27700) %>%
    rename(geometry=x) %>%
    mutate(siteid=filter(spoke.df, side=="start")$siteid)
  sf.quad <- sf$siteid %>%
    map_dfr(~st_split(filter(sf, siteid==.x), filter(spoke.lines, siteid==.x)) %>%
              st_collection_extract() %>%
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
  combos <- crossing(indices, var=as.list(variable))
  
  quosures <- map2(combos$indices, combos$var,
                   ~quo(lag(!!.y, !!.x)) ) %>%
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
  y.df %>%
    rowwise() %>%
    mutate(tl=tl_i$tl[max(which(y==tl_i$abbr & {{N}} >= tl_i$min_ge))],
           alert=tl_i$A[max(which(y==tl_i$abbr & {{N}} >= tl_i$min_ge))]) %>%
    ungroup# %>%
    # mutate(alert=c("0_none", "1_warn", "2_alert")[alert+1])
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
calc_y_features <- function(yRaw.df, y_i, tl_i, dist.df=NULL) {
  y.ls <- yRaw.df %>% 
    group_by(siteid, date) %>%
    slice_head(n=1) %>%
    ungroup %>%
    pivot_longer(any_of(y_i$abbr), names_to="y", values_to="N") %>%
    filter(!is.na(N)) %>%
    mutate(lnN=log1p(N)) %>%
    get_trafficLights(N, tl_i) %>%
    arrange(y, siteid, date) %>%
    group_by(y, siteid) %>%
    get_lags(lnN, alert, date, n=2) %>%
    ungroup %>%
    mutate(year=year(date),
           lnNWt1=lnN1/log1p(as.numeric(date-date1)), 
           lnNWt2=lnN2/log1p(as.numeric(date-date2)),
           lnNAvg1=0, lnNAvg2=0, prAlertAvg1=0, prAlertAvg2=0, 
           lnNAvgPrevYr=0, prAlertAvgPrevYr=0) %>%
    group_by(y) %>%
    group_split()
  y_new <- vector("list", length(y.ls))
  
  for(i in seq_along(y.ls)) {
    y.df_i <- y.ls[[i]] %>% select(siteid, date, lnN1, lnN2, alert1, alert2)
    yYr_i <- y.ls[[i]] %>% 
      group_by(siteid, year) %>%
      summarise(lnN_yr=mean(lnN, na.rm=T),
                prAlert_yr=mean(alert=="A1", na.rm=T)) %>%
      ungroup %>%
      mutate(year_matchObs=year+1)
    y_new[[i]] <- y.ls[[i]] %>%
      left_join(yYr_i %>% select(year_matchObs, siteid, lnN_yr, prAlert_yr), 
                by=c("year"="year_matchObs", "siteid"="siteid")) %>%
      rename(lnNPrevYr=lnN_yr, prAlertPrevYr=prAlert_yr)
    
    for(j in 1:nrow(y.ls[[i]])) {
      site_j <- y.df_i$siteid[j]
      date_j <- y.df_i$date[j]
      yr_j <- year(date_j)
      wk.df <- y.df_i %>%
        filter(siteid %in% dist.df$dest_c[dist.df$origin==site_j][[1]] &
                 date <= date_j & 
                 date > date_j-7) 
      yr.df <- yYr_i %>%
        filter(siteid %in% dist.df$dest_c[dist.df$origin==site_j][[1]] &
                 year == yr_j - 1)
      yr.site_j <- yYr_i %>% filter(siteid==site_j & year==yr_j-1)
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
  hab_ids <- site_tox.sf %>%
    select(siteid, geom) %>%
    mutate(hab_id=st_intersects(., site_hab.sf)) %>%
    st_drop_geometry()
  
  # match hab observations to toxin observations
  hab.df <- hab.df %>%
    mutate(prA=as.numeric(alert=="A1")) %>%
    select(siteid, date, y, lnN, prA) %>% 
    rename(hab_id=siteid, lnNAvg=lnN) %>%
    pivot_wider(names_from=y, values_from=c(lnNAvg, prA), names_glue="{y}{.value}")
  hab_y_names <- names(hab.df)[-(1:2)]
  habSums <- vector("list", nrow(tox.obs))
  for(i in 1:nrow(tox.obs)) {
    hab_sites <- filter(hab_ids, siteid==tox.obs$siteid[i])$hab_id[[1]]
    habSums[[i]] <- hab.df %>% 
      filter(hab_id %in% hab_sites & 
               date <= tox.obs$date[i] - 7*1 &
               date >= tox.obs$date[i] - 7*5) %>%
      summarise(across(all_of(hab_y_names), ~mean(.x, na.rm=T)))
  }
  out.df <- tox.obs %>% bind_cols(do.call('rbind', habSums))
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
  
  d.ls$compiled <- d.ls$site %>% select(-sin) %>%
    right_join(d.ls$obs, by="siteid", multiple="all") %>%
    left_join(d.ls$cmems.pt, by=c("cmems_id", "date")) %>%
    left_join(d.ls$cmems.buf %>%
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
              by=c("siteid", "date")) %>%
    mutate(wrf_id=if_else(date < "2019-04-01", wrf_id.1, wrf_id.2)) %>%
    select(-wrf_id.1, -wrf_id.2, -version) %>%
    left_join(d.ls$wrf.pt %>% select(-version), by=c("wrf_id", "date")) %>%
    left_join(d.ls$wrf.buf %>%
                pivot_wider(names_from="quadrant", values_from=-(1:3), names_sep="Dir"),
              by=c("siteid", "date")) %>%
    mutate(year=year(date),
           yday=yday(date))
  if(target=="tox") {
    d.ls$compiled <- d.ls$compiled %>%
      left_join(d.ls$habAvg %>% select(-date, -siteid), by="obsid")
  }
  d.ls$compiled <- d.ls$compiled %>% 
    na.omit
  
  return(d.ls)
}






# Model preparation -------------------------------------------------------


#' Create recipe and prepare using training data
#'
#' @param train.df 
#' @param response 
#' @param dimReduce 
#'
#' @return
#' @export
#'
#' @examples
prep_recipe <- function(train.df, response, covsExclude="NA", dimReduce=FALSE) {
  respExclude <- grep(response, c("lnN", "tl", "alert"), value=T, invert=T)
  pred_vars <- names(train.df)
  include_interactions <- !grepl("X", covsExclude)
  rec <- recipe(train.df) %>%
    step_mutate(prevAlert=alert1, role="ID") %>%
    update_role(all_of(pred_vars), new_role="predictor") %>%
    update_role(all_of(response), new_role="outcome") %>%
    update_role(obsid, y, date, siteid, year, new_role="ID") %>%
    update_role(lon, lat, new_role="RE") %>%
    step_select(-any_of(respExclude)) %>%
    step_dummy(all_factor_predictors()) %>%
    step_logit(starts_with("prAlert"), offset=0.01) %>%
    step_logit(ends_with("A1"), offset=0.01) %>%
    step_harmonic(yday, frequency=1, cycle_size=365, keep_original_cols=T) %>%
    update_role(yday, new_role="RE") %>%
    step_rename(ydayCos=yday_cos_1, ydaySin=yday_sin_1) %>%
    step_mutate_at(lon, lat, fn=list(z=~.)) %>%
    step_interact(term=~ydaySin:ydayCos, sep="X")
  if(include_interactions) {
    rec <- rec %>%
      step_interact(terms=~UWk:fetch:matches("Dir[EW]"), sep="X") %>%
      step_interact(terms=~VWk:fetch:matches("Dir[NS]"), sep="X")
  }
  rec <- rec %>%
    step_interact(terms=~lon_z:lat_z, sep="X") %>%
    step_YeoJohnson(all_predictors()) %>%
    step_normalize(all_predictors()) %>%
    step_rename_at(contains("_"), fn=~str_remove_all(.x, "_")) %>%
    step_select(-matches(covsExclude))
  if(dimReduce) {
    rec <- rec %>%
      step_pca(all_predictors(), threshold=0.95)
  }
  rec %>%
    prep(training=train.df)
}






#' Filter list of covariates following recipe thinning
#'
#' @param all_covs 
#' @param data.y 
#' @param covsExclude 
#'
#' @return
#' @export
#'
#' @examples
filter_corr_covs <- function(all_covs, data.y, covsExclude="NA") {
  uncorr_covs <- unique(unlist(map(data.y, ~unlist(map(.x, names)))))
  if(any(grepl("^PC", uncorr_covs))) {
    covs_incl <- list(main=grep("^PC", uncorr_covs, value=T),
                      interact=NULL,
                      nonHB=NULL)
  } else {
    uncorr_covs <- unique(c(uncorr_covs, "ydayCos", "ydaySin", "lon", "lat"))
    covs_incl <- map(all_covs, ~.x[.x %in% uncorr_covs])
  }
  
  covs_incl <- map(covs_incl, ~grep(covsExclude, .x, value=T, invert=T))
  return(covs_incl)
}





#' Create formulas for Hierarchical Bayesian models
#'
#' @param resp 
#' @param covs 
#' @param sTerms 
#' @param splinesInt 
#' @param splinesCovs 
#'
#' @return
#' @export
#'
#' @examples
make_HB_formula <- function(resp, covs, sTerms=NULL, 
                            splinesInt="both", splinesCovs="time") {
  library(tidyverse); library(brms); library(glue)
  
  splines_int <- switch(splinesInt,
                        "time"="s(yday, bs=('cc'))",
                        "space"="s(lon, lat)",
                        "both"="t2(yday, lon, lat, bs=c('cc','ts'), d=c(1,2))")
  splines_cov <- switch(splinesCovs,
                        "time"="s(yday, bs=('cc'))",
                        "space"="s(lon, lat)",
                        "both"="t2(yday, lon, lat, bs=c('cc','ts'), d=c(1,2))")
  
  if(is.null(sTerms)) {
    if(resp=="lnN") {
      form <- bf(glue("{resp} ~ 1 + {paste(covs, collapse='+')}", 
                      "+ {splines_int}",
                      "+ (1 + {paste(covs, collapse='+')} | siteid)"),
                 glue("hu ~ 1 + {paste(covs, collapse='+')}", 
                      "+ {splines_int}",
                      "+ (1 + {paste(covs, collapse='+')} | siteid)"))  
    } else {
      form <- bf(glue("{resp} ~ 1 + {paste(covs, collapse='+')}", 
                      "+ {splines_int}",
                      "+ (1 + {paste(covs, collapse='+')} | siteid)"))
    }
  } else {
    flist <- c(
      map(sTerms$b, ~as.formula(glue("{.x} ~ {splines_cov} + (1|siteid)"))),
      map(sTerms$p, ~as.formula(glue("{.x} ~ 1 + (1|siteid)"))),
      map(1, ~glue("bIntercept ~ 1 + {splines_int} + (1|siteid)"))
    )
    if(resp=="lnN") {
      sTerms$pHu <- paste0(sTerms$p, "Hu")
      sTerms$bHu <- paste0(sTerms$b, "Hu")
      flist <- c(
        nlf(glue("hu ~ 1*bInterceptHu",
                 "+ {paste(sTerms$pHu, sTerms$bHu, covs, sep='*', collapse='+')}")),
        flist,
        map(sTerms$bHu, ~as.formula(glue("{.x} ~ {splines_cov} + (1|siteid)"))),
        map(sTerms$pHu, ~as.formula(glue("{.x} ~ 1 + (1|siteid)"))),
        map(1, ~glue("bInterceptHu ~ 1 + {splines_int} + (1|siteid)"))
      )
      form <- bf(glue("{resp} ~ 1*bIntercept",
                      "+ {paste(sTerms$p, sTerms$b, covs, sep='*', collapse='+')}"),
                 flist=flist, nl=T) 
    } else {
      form <- bf(glue("{resp} ~ 1*bIntercept",
                      "+ {paste(sTerms$p, sTerms$b, covs, sep='*', collapse='+')}"),
                 flist=flist, nl=T)
    }
  }
  return(form)
}







#' Create priors for each Hierarchical Bayesian model
#'
#' @param prior_i 
#' @param mod 
#' @param resp 
#' @param covs 
#'
#' @return
#' @export
#'
#' @examples
make_HB_priors <- function(prior_i, mod, resp, covs, PCA=F) {
  library(tidyverse); library(brms)
  if(mod=="HBL") {
    p <- c(prior(normal(0, 1), class="Intercept"),
           prior(normal(0, 0.1), class="sd"))
    if(resp=="tl") {
      p <- c(p,
             prior_string(glue("horseshoe({prior_i$hs1}, par_ratio={prior_i$hs2})"), class="b"))
    } else {
      p <- c(p,
             prior_string(glue("R2D2({prior_i$r1}, {prior_i$r2})"), class="b"))
    }
    if(resp=="lnN") {
      p <- c(p,
             prior_string(glue("R2D2({prior_i$r1}, {prior_i$r2})"), class="b", dpar="hu"))
    }
  }
  if(mod=="HBN") {
    if(PCA) {
      terms <- list(int="bIntercept",
                    cov=covs)
    } else {
      terms <- list(int="bIntercept",
                    cov=c(covs$main, covs$interact))
    }
    
    if(resp=="lnN") {
      terms <- map(terms, ~c(.x, paste0(.x, "Hu")))
    }
    p <- c(
      map(terms$int, 
          ~c(prior_string("normal(0,1)", class="b", nlpar=.x),
             prior_string("normal(0,.5)", class="sd", nlpar=.x, lb=0),
             prior_string("student_t(3,0,2.5)", class="sds", nlpar=.x, lb=0))) %>%
        do.call('c', .),
      map(terms$cov, 
          ~c(prior_string(glue("beta({prior_i$b},1)"), 
                          nlpar=paste0("p", .x), lb=0, ub=1),
             prior_string("normal(0,1)", class="b", 
                          nlpar=paste0("b", .x)),
             prior_string("normal(0,.5)", class="sd", 
                          nlpar=paste0("b", .x), lb=0),
             prior_string("normal(0,.5)", class="sd", 
                          nlpar=paste0("p", .x), lb=0),
             prior_string(glue("double_exponential(0,{prior_i$de})"), class="sds", 
                          nlpar=paste0("b", .x), lb=0))) %>%
        do.call('c', .)
    )
  }
  return(p)
}





#' Identify indexes for ML training
#'
#' @param data.df 
#'
#' @return
#' @export
#'
#' @examples
createFoldsByYear <- function(data.df) {
  folds_out <- data.df %>% mutate(rowNum=row_number()) %>% 
    group_by(year) %>% group_split() %>% map(~.x$rowNum)
  folds_out <- folds_out[-1]
  folds_in <- map(folds_out, ~(1:nrow(data.df))[-.x])
  return(list(i.in=folds_in, i.out=folds_out))
}







#' Wrapper to fit a model
#'
#' @param mod 
#' @param resp 
#' @param form.ls 
#' @param d.ls 
#' @param opts 
#' @param tunes 
#' @param out.dir 
#' @param y 
#' @param suffix 
#'
#' @return
#' @export
#'
#' @examples
fit_model <- function(mod, resp, form.ls, d.ls, opts, tunes, out.dir, y, suffix=NULL) {
  library(glue); library(tidymodels)
  dir.create(glue("{out.dir}/meta/"), showWarnings=F, recursive=T)
  PCA_run <- all(!is.null(suffix), grepl("PCA", suffix))
  # Fit ML models
  if(mod %in% c("Ridge", "ENet", "RF", "NN", "MARS", "Boost")) {
    fit_ID <- glue("{y}_{resp}_{mod}{ifelse(is.null(suffix),'',suffix)}")
    mod.prefix <- ifelse(PCA_run, "PCA.", "")
    ML_form <- ifelse(PCA_run, "ML_PCA", "ML")
    ML_spec <- switch(mod,
                      Ridge=logistic_reg(penalty=tune(), 
                                         mixture=0) %>%
                        set_engine("glmnet") %>% set_mode("classification"),
                      ENet=logistic_reg(penalty=tune(), 
                                        mixture=tune()) %>%
                        set_engine("glmnet") %>% set_mode("classification"),
                      RF=rand_forest(trees=tune(), 
                                     min_n=tune()) %>%
                        set_engine("randomForest") %>% set_mode("classification"),
                      NN=mlp(hidden_units=tune(), 
                             penalty=tune(), 
                             epochs=tune()) %>%
                        set_engine("nnet", maxNWts=1e4) %>% set_mode("classification"),
                      MARS=mars(num_terms=tune(), 
                                prod_degree=tune()) %>%
                        set_engine("earth") %>% set_mode("classification"),
                      Boost=boost_tree(trees=tune(),
                                       tree_depth=tune(),
                                       min_n=tune(),
                                       learn_rate=tune(),
                                       loss_reduction=tune()) %>%
                        set_engine("xgboost") %>% set_mode("classification")
    )
    avg_prec2 <- metric_tweak("avg_prec2", average_precision, event_level="second")
    pr_auc2 <- metric_tweak("pr_auc2", pr_auc, event_level="second")
    roc_auc2 <- metric_tweak("roc_auc2", roc_auc, event_level="second")
    wf <- workflow() %>%
      add_model(ML_spec) %>%
      add_formula(form.ls[[resp]][[ML_form]])
    set.seed(1003)
    out_tune <- wf %>%
      tune_grid(resamples=opts, 
                grid=grid_latin_hypercube(extract_parameter_set_dials(ML_spec), 
                                          size=tunes[[mod]]),
                metrics=metric_set(roc_auc2, pr_auc2, avg_prec2), 
                control=control_grid(save_pred=T))
    saveRDS(out_tune %>% butcher, glue("{out.dir}/meta/{fit_ID}_tune.rds"))
    best <- select_best(out_tune, "avg_prec2")
    out_tune %>% 
      collect_predictions() %>%
      filter(.config==best$.config) %>%
      arrange(.row) %>%
      mutate(obsid=d.ls[[resp]]$obsid,
             y=y) %>%
      select(y, obsid, .pred_A1) %>%
      rename_with(~glue("{mod.prefix}{mod}_{resp}_A1"), .cols=".pred_A1") %>%
      saveRDS(glue("{out.dir}/cv/{fit_ID}_CV.rds"))
    out <- wf %>%
      finalize_workflow(best) %>%
      fit(data=d.ls[[resp]]) %>%
      butcher()
  }
  
  # Fit Hierarchical Bayesian models
  if(mod %in% c("HBL", "HBN")) {
    library(brms)
    fit_ID <- glue("{y}_{resp}_{mod}{opts$prior_i}{ifelse(is.null(suffix),'',suffix)}")
    HB.family <- switch(resp, 
                        lnN=hurdle_lognormal,
                        tl=cumulative,
                        alert=bernoulli)
    HB_form <- ifelse(PCA_run, paste0(mod, "_PCA"), paste0(mod))
    HB_form_dummy <- ifelse(PCA_run, "HB_vars_PCA", "HB_vars")
    wf <- workflow() %>%
      add_model(bayesian(mode="classification", engine="brms", 
                         formula.override=bayesian_formula(form.ls[[resp]][[HB_form]]),
                         family=HB.family, 
                         prior=tunes[[resp]][[HB_form]],
                         init=0, 
                         iter=opts$iter,
                         warmup=opts$warmup,
                         control=opts$ctrl,
                         chains=opts$chains,
                         cores=opts$cores,
                         save_model=glue("{out.dir}/meta/{fit_ID}.stan")),
                formula=form.ls[[resp]]$HB_vars) %>%
      add_recipe(recipe(d.ls[[resp]], formula=form.ls[[resp]][[HB_form_dummy]]))
    out <- wf %>%
      fit(data=d.ls[[resp]]) %>%
      axe_env_bayesian()
  }
  saveRDS(out, glue("{out.dir}/{fit_ID}.rds"))
  cat("Saved ", y, "_", resp, "_", mod, " as ", out.dir, "*", suffix, "\n", sep="")
}









#' Run cross-validation for Bayesian models
#'
#' @param mod 
#' @param folds 
#' @param cv.dir 
#' @param y 
#' @param y_i.i 
#' @param r 
#' @param form.ls 
#' @param HB.i 
#' @param priors 
#'
#' @return
#' @export
#'
#' @examples
run_Bayes_CV <- function(mod, folds, cv.dir, y, y_i.i, r, form.ls, HB.i, priors, PCA=F) {
  for(f in 1:nrow(folds)) {
    f_ <- paste0("_f", str_pad(f, 2, side="left", pad="0"))
    f_ <- ifelse(PCA, paste0("_PCA", f_), f_)
    d.cv <- list(train=list(alert=training(folds$splits[[f]])),
                 test=list(alert=testing(folds$splits[[f]])))
    if(!file.exists(glue("{cv.dir}/{y}_{r}_{mod}_CV{f_}.rds"))) {
      fit_model(mod, r, form.ls, d.cv$train, HB.i, priors, cv.dir, y, f_)
    }
    if(file.exists(glue("{cv.dir}/{y}_{r}_{mod}1{f_}.rds"))) {
      summarise_predictions(d.cv$test, d.cv$test, r, cv.dir, y_i.i, f_) %>%
        saveRDS(glue("{cv.dir}/{y}_{r}_{mod}_CV{f_}.rds"))
      file.remove(glue("{cv.dir}/{y}_{r}_{mod}1{f_}.rds"))
    }
  }
}






# Model output ------------------------------------------------------------





#' Summarise and aggregate predictions for a set of models
#'
#' @param d.y 
#' @param set 
#' @param resp 
#' @param fit.dir 
#' @param y_i.i 
#' @param suffix 
#'
#' @return
#' @export
#'
#' @examples
summarise_predictions <- function(d.y, dPCA.y, resp, fit.dir, y_i.i, suffix=NULL) {
  library(tidyverse); library(glue); library(tidymodels)
  fits.f <- dirf(fit.dir, glue("{y_i.i$abbr[1]}_{resp}.*{ifelse(is.null(suffix),'',suffix)}"))
  names(fits.f) <- str_split_fixed(str_split_fixed(fits.f, glue("{resp}_"), 2)[,2], "_|\\.", 2)[,1]
  fits <- map(fits.f, readRDS)
  
  fits.dPCA <- grep("PCA", fits.f)
  fits.d <- grep("PCA", fits.f, invert=T)
  
  if(length(fits.dPCA) > 0) {
    preds.dPCA <- imap_dfc(fits[fits.dPCA], ~get_predictions(.x, .y, resp, dPCA.y, y_i.i)) %>%
      rename_with(~glue("PCA.{.x}"))
  } else {
    preds.dPCA <- NULL
  }
  
  if(length(fits.d) > 0) {
    preds.d <- imap_dfc(fits[fits.d], ~get_predictions(.x, .y, resp, d.y, y_i.i))
  } else {
    preds.d <- NULL
  }
  
  out.df <- d.y[[resp]] %>%
    select(y, obsid, siteid, date, {{resp}}) 
  if(!is.null(preds.d)) {
    out.df <- out.df %>% bind_cols(preds.d)
  }
  if(!is.null(preds.dPCA)) {
    out.df <- out.df %>% bind_cols(preds.dPCA)
  }
  return(out.df)
}





#' Generate predictions for a model
#'
#' @param fit 
#' @param mod 
#' @param resp 
#' @param d.df 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
get_predictions <- function(fit, mod, resp, d.df, y_i.i) {
  library(tidyverse); library(glue); library(tidymodels); library(brms)
  
  if(grepl("HB", mod)) {
    pred <- parsnip::extract_fit_engine(fit) %>%
      posterior_epred(d.df[[resp]], allow_new_levels=T) %>%
      summarise_post_preds(., resp, y_i.i)
  } else {
    pred_type <- ifelse(resp=="lnN", "raw", "prob")
    preds <- predict(fit, d.df[[resp]], pred_type)
    pred <- summarise_ML_preds(preds, resp, y_i.i)
  }
  pred.df <- as_tibble(pred) %>% 
    rename_with(~glue("{mod}_{resp}_{.x}"))
  return(pred.df)
}




#' Calculate mean predictions for ML models
#'
#' @param preds 
#' @param resp 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
summarise_ML_preds <- function(preds, resp, y_i.i) {
  if(resp=="alert") {
    pred <- cbind(A1=preds$.pred_A1)
  }
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1)) + 1
    pred <- cbind(preds, 
                  rowSums(preds[,thresh:4]))
    colnames(pred) <- c(paste0("TL", 0:3), "A1")
  }
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    pred <- cbind(lnN=preds,
                  A1=preds>=thresh)
  }
  return(pred)
}




#' Calculate mean predictions for HB models
#'
#' @param post 
#' @param resp 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
summarise_post_preds <- function(post, resp, y_i.i) {
  if(resp=="alert") {
    pred <- cbind(A1=colMeans(post))
  }
  
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1)) + 1
    pred <- cbind(apply(post, 2:3, mean),
                  colMeans(apply(post[,,(thresh):4], 1:2, sum)))
    colnames(pred) <- c(paste0("TL", 0:3), "A1")
  }
  
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    pred <- cbind(lnN=colMeans(post),
                  A1=colMeans(post>=thresh))
  }
  return(pred)
} 









#' Merge summarised predictions across all models into single dataframes
#'
#' @param files vector of full file path
#'
#' @return
#' @export
#'
#' @examples
merge_pred_dfs <- function(files, CV=NULL) {
  f.df <- tibble(f=files, 
                 covSet=str_split(files, "/") %>% 
                   map_chr(~grep("^[1-9]", .x, value=T) %>% str_sub(1, 1)))
  if(is.null(CV)) {
    map(1:nrow(f.df), 
        ~readRDS(f.df$f[.x]) %>% 
          lapply(., function(x) {x %>% mutate(covSet=paste0("d", f.df$covSet[.x], "."))})) %>%
      list_transpose() %>%
      map_depth(2, ~.x %>% 
                  pivot_longer(ends_with("_A1"), names_to="model", values_to="prA1") %>%
                  mutate(model=paste0(covSet, model)) %>%
                  select(-covSet) %>%
                  pivot_wider(names_from="model", values_from="prA1")) %>%
      map(~reduce(.x, full_join)) 
  } else if(CV=="HB") {
    map_dfr(1:nrow(f.df), 
            ~readRDS(f.df$f[.x]) %>% mutate(covSet=paste0("d", f.df$covSet[.x], "."))) %>%
      pivot_longer(ends_with("A1"), names_to="model", values_to="prA1") %>%
      mutate(model=paste0(covSet, model)) %>%
      select(-covSet) %>%
      pivot_wider(names_from="model", values_from="prA1")
  } else if(CV=="ML") {
    map(1:nrow(f.df), 
        ~readRDS(f.df$f[.x]) %>% 
          mutate(covSet=paste0("d", f.df$covSet[.x], ".")) %>%
          pivot_longer(ends_with("A1"), names_to="model", values_to="prA1") %>%
          mutate(model=paste0(covSet, model)) %>%
          select(-covSet) %>%
          pivot_wider(names_from="model", values_from="prA1")) %>%
      reduce(full_join)
  } else {
    stop("CV must be 'HB', 'ML', or NULL")
  }
}








#' Calculate model weights based on mean log loss
#'
#' @param cv.df 
#' @param resp 
#' @param wt.penalty 
#'
#' @return
#' @export
#'
#' @examples
calc_LL_wts <- function(cv.df, resp, wt.penalty=2) {
  library(yardstick)
  if(resp=="alert") {
    wt.df <- cv.df %>%
      pivot_longer(ends_with("_A1"), names_to="model", values_to="pr") %>%
      group_by(model) %>%
      average_precision(pr, truth=alert, event_level="second") %>%
      mutate(.estimate=log(.estimate))
  }
  if(resp=="tl") {
    wt.df <- cv.df %>%
      pivot_longer(contains("_tl_TL"), names_to="ID", values_to="pr") %>%
      mutate(model=str_split_fixed(ID, "_TL", 2)[,1],
             predCat=str_split_fixed(ID, "_tl_", 2)[,2]) %>%
      select(-ID) %>%
      pivot_wider(names_from=predCat, values_from=pr) %>%
      group_by(model) %>%
      mn_log_loss(TL0, TL1, TL2, TL3, truth=tl)
  }
  if(resp=="lnN") {
    wt.df <- cv.df %>%
      pivot_longer(ends_with("_lnN_lnN"), names_to="model", values_to="pr") %>%
      group_by(model) %>%
      rmse(truth=lnN, estimate=pr)
  }
  return(wt.df %>%
           ungroup %>%
           mutate(wt=(1/.estimate^wt.penalty)/sum(1/.estimate^wt.penalty)))
}





#' Calculate ensemble model predictions
#'
#' @param out.ls 
#' @param wt.ls 
#' @param resp 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
calc_ensemble <- function(out.ls, wt.ls, resp, y_i.i, method="wtmean", out.path=NULL, opt=NULL) {
  library(tidyverse); library(tidymodels)
  if(grepl("RF|GLM|HB", method) & resp != "alert") {
    stop("Models only implemented for alert")
  }
  if(resp=="alert") {
    if(method=="wtmean") {
      out <- left_join(
        out.ls[[resp]] %>% 
          pivot_longer(ends_with("_A1"), names_to="model", values_to="pr"),
        wt.ls[[resp]]
      ) %>%
        mutate(pr_logit=brms::logit_scaled(pr, lb=-0.01, ub=1.01)) %>%
        group_by(obsid) %>%
        summarise(ens_alert_A1=sum(pr*wt, na.rm=T),
                  ensLogitMn_alert_A1=brms::inv_logit_scaled(
                    sum(pr_logit*wt, na.rm=T))) %>%
        ungroup
    } else if(grepl("[GLM|RF|HB]_fit", method)) {
      avg_prec2 <- metric_tweak("avg_prec2", average_precision, event_level="second")
      folds <- vfold_cv(wt.ls[[resp]], strata="alert")
      ens_rec <- recipe(alert~., data=wt.ls[[resp]]) %>%
        update_role(y, obsid, siteid, date, new_role="ID") %>%
        step_logit(ends_with("_A1"), offset=0.01) %>%
        step_normalize(all_predictors())
      ens_rec2 <- recipe(alert~., data=wt.ls[[resp]]) %>%
        update_role(y, obsid, siteid, date, new_role="ID") 
      
      if(method=="GLM_fit") {
        size <- ifelse(is.null(opt), 1e3, opt)
        GLM_spec <- logistic_reg(penalty=tune(), mixture=0) %>%
          set_engine("glmnet", lower.limits=0) %>% set_mode("classification")
        GLM_wf <- workflow() %>%
          add_model(GLM_spec) %>%
          add_recipe(ens_rec)
        GLM_tune <- GLM_wf %>%
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(GLM_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        GLM_out <- GLM_wf %>%
          finalize_workflow(select_best(GLM_tune, "avg_prec2")) %>%
          fit(wt.ls[[resp]]) %>%
          butcher
        saveRDS(GLM_out, glue("{out.path}/{y_i.i$abbr}_EnsGLM.rds"))
      }
      
      if(method=="RF_fit") {
        size <- ifelse(is.null(opt), 1e2, opt)
        RF_spec <- rand_forest(trees=tune(), 
                               min_n=tune()) %>%
          set_engine("randomForest") %>% set_mode("classification")
        RF_wf <- workflow() %>%
          add_model(RF_spec) %>%
          add_recipe(ens_rec)
        RF_tune <- RF_wf %>%
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(RF_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        RF_out <- RF_wf %>%
          finalize_workflow(select_best(RF_tune, "avg_prec2")) %>%
          fit(wt.ls[[resp]]) %>%
          butcher
        saveRDS(RF_out, glue("{out.path}/{y_i.i$abbr}_EnsRF.rds"))
      }
      
      if(method=="HB_fit") {
        mn <- ifelse(is.null(opt), 0.5, opt)
        library(brms)
        wf1 <- workflow() %>%
          add_model(bayesian(mode="classification", engine="brms", 
                             family=bernoulli, 
                             prior=prior_string(glue("R2D2({mn})"), class="b"),
                             init=0, 
                             iter=1500,
                             warmup=1000,
                             chains=3,
                             cores=3,
                             save_model=glue("{out.path}/{y_i.i$abbr}_EnsHB1.stan")),
                    formula=alert~.) %>%
          add_recipe(ens_rec)
        HB_out1 <- wf1 %>%
          fit(data=wt.ls[[resp]]) 
        saveRDS(HB_out1, glue("{out.path}/{y_i.i$abbr}_EnsHB1.rds"))
        
        wf2 <- workflow() %>%
          add_model(bayesian(mode="classification", engine="brms", 
                             family=bernoulli, 
                             prior=prior(exponential(2), class=b, lb=0),
                             init=0, 
                             iter=1500,
                             warmup=1000,
                             chains=3,
                             cores=3,
                             save_model=glue("{out.path}/{y_i.i$abbr}_EnsHB2.stan")),
                    formula=alert~.) %>%
          add_recipe(ens_rec)
        HB_out2 <- wf2 %>%
          fit(data=wt.ls[[resp]]) 
        saveRDS(HB_out2, glue("{out.path}/{y_i.i$abbr}_EnsHB2.rds"))
        
        wf3 <- workflow() %>%
          add_model(bayesian(mode="classification", engine="brms", 
                             family=bernoulli, 
                             prior=prior_string(glue("R2D2({mn})"), class="b"),
                             init=0, 
                             iter=1500,
                             warmup=1000,
                             chains=3,
                             cores=3,
                             save_model=glue("{out.path}/{y_i.i$abbr}_EnsHB3.stan")),
                    formula=alert~.) %>%
          add_recipe(ens_rec2)
        HB_out3 <- wf3 %>%
          fit(data=wt.ls[[resp]]) 
        saveRDS(HB_out3, glue("{out.path}/{y_i.i$abbr}_EnsHB3.rds"))
        
        wf4 <- workflow() %>%
          add_model(bayesian(mode="classification", engine="brms", 
                             family=bernoulli, 
                             prior=prior(exponential(2), class=b, lb=0),
                             init=0, 
                             iter=1500,
                             warmup=1000,
                             chains=3,
                             cores=3,
                             save_model=glue("{out.path}/{y_i.i$abbr}_EnsHB4.stan")),
                    formula=alert~.) %>%
          add_recipe(ens_rec2)
        HB_out4 <- wf4 %>%
          fit(data=wt.ls[[resp]]) 
        saveRDS(HB_out4, glue("{out.path}/{y_i.i$abbr}_EnsHB4.rds"))
      }
    }
    if(grepl("GLM", method)) {
      GLM_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsGLM.rds"))
      out <- out.ls[[resp]] %>%
        mutate(ensGLM_alert_A1=predict(GLM_out, new_data=., type="prob")[[2]]) %>%
        select(obsid, ensGLM_alert_A1) 
    }
    if(grepl("RF", method)) {
      RF_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsRF.rds"))
      out <- out.ls[[resp]] %>%
        mutate(ensRF_alert_A1=predict(RF_out, new_data=., type="prob")[[2]]) %>%
        select(obsid, ensRF_alert_A1) 
    }
    if(grepl("HB", method)) {
      HB_out1 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsHB1.rds"))
      HB_out2 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsHB2.rds"))
      HB_out3 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsHB3.rds"))
      HB_out4 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsHB4.rds"))
      pred1 <- parsnip::extract_fit_engine(HB_out1) %>%
        posterior_epred(out.ls[[resp]], allow_new_levels=T) %>%
        summarise_post_preds(., resp, y_i.i)
      pred2 <- parsnip::extract_fit_engine(HB_out2) %>%
        posterior_epred(out.ls[[resp]], allow_new_levels=T) %>%
        summarise_post_preds(., resp, y_i.i)
      pred3 <- parsnip::extract_fit_engine(HB_out3) %>%
        posterior_epred(out.ls[[resp]], allow_new_levels=T) %>%
        summarise_post_preds(., resp, y_i.i)
      pred4 <- parsnip::extract_fit_engine(HB_out4) %>%
        posterior_epred(out.ls[[resp]], allow_new_levels=T) %>%
        summarise_post_preds(., resp, y_i.i)
      out <- out.ls[[resp]] %>%
        mutate(ensHB1_alert_A1=pred1[,1],
               ensHB2_alert_A1=pred2[,1],
               ensHB3_alert_A1=pred3[,1],
               ensHB4_alert_A1=pred4[,1]) %>%
        select(obsid, starts_with("ensHB")) 
    }
  }
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1))
    alert_cols <- paste0("ens_tl_TL", thresh:3)
    out <- left_join(
      out.ls[[resp]] %>%
        select(-ends_with("_A1")) %>%
        pivot_longer(contains("_tl_TL"), names_to="ID", values_to="pr") %>%
        mutate(model=str_split_fixed(ID, "_TL", 2)[,1],
               predCat=str_split_fixed(ID, "_tl_", 2)[,2]) %>%
        select(-ID),
      wt.ls[[resp]]) %>%
      group_by(obsid, predCat) %>%
      summarise(ens_tl=sum(pr*wt, na.rm=T)) %>%
      ungroup %>%
      pivot_wider(names_from=predCat, values_from=ens_tl, names_prefix="ens_tl_") %>%
      mutate(ens_tl_A1=rowSums(pick(all_of(alert_cols))))
  }
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    out <- left_join(
      out.ls[[resp]] %>%
        select(-ends_with("_A1")) %>%
        pivot_longer(ends_with("_lnN_lnN"), names_to="model", values_to="pr"),
      wt.ls[[resp]]) %>%
      group_by(obsid) %>%
      summarise(ens_lnN=sum(pr*wt, na.rm=T),
                ens_lnN_A1=sum((pr >= thresh)*wt, na.rm=T)) %>%
      ungroup
  }
  return(left_join(out.ls[[resp]], out))
}






#' Identify rows in y within thresh days of x 
#'
#' @param x 
#' @param y 
#' @param thresh 
#'
#' @return
#' @export
#'
#' @examples
get_rows_by_yday <- function(x, y, thresh) {
  xy_diff <- x - y
  comp.mx <- cbind(abs(xy_diff),
                   abs(xy_diff - 365),
                   abs(xy_diff + 365))
  ids <- which(apply(comp.mx, 1, function(j) any(j <= thresh)))
  return(ids)
}






#' Calculate null model predictions
#'
#' @param obs.ls 
#' @param resp 
#'
#' @return
#' @export
#'
#' @examples
calc_null <- function(obs.ls, resp) {
  library(tidyverse)
  obs.df <- obs.ls[[resp]] %>%
    mutate(yday=yday(date)) %>%
    select(yday, {{resp}})
  if(resp=="alert") {
    temp.df <- tibble(yday=1:365,
                      null4wk_alert_A1=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_alert_A1[i] <- mean(obs.df$alert[obs.ids] == "A1")
    }
    out.df <- obs.ls[[resp]] %>% 
      mutate(yday=yday(date),
             nullGrand_alert_A1=mean(alert=="A1"))
  }
  if(resp=="tl") {
    temp.df <- tibble(yday=1:365,
                      null4wk_tl_TL0=NA_real_,
                      null4wk_tl_TL1=NA_real_,
                      null4wk_tl_TL2=NA_real_,
                      null4wk_tl_TL3=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_tl_TL0[i] <- mean(obs.df$tl[obs.ids] == "TL0")
      temp.df$null4wk_tl_TL1[i] <- mean(obs.df$tl[obs.ids] == "TL1")
      temp.df$null4wk_tl_TL2[i] <- mean(obs.df$tl[obs.ids] == "TL2")
      temp.df$null4wk_tl_TL3[i] <- mean(obs.df$tl[obs.ids] == "TL3")
    }
    out.df <- obs.ls[[resp]] %>% 
      mutate(yday=yday(date),
             nullGrand_tl_TL0=mean(tl=="TL0"),
             nullGrand_tl_TL1=mean(tl=="TL1"),
             nullGrand_tl_TL2=mean(tl=="TL2"),
             nullGrand_tl_TL3=mean(tl=="TL3"))
  }
  if(resp=="lnN") {
    temp.df <- tibble(yday=1:365,
                      null4wk_lnN_lnN=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_lnN_lnN[i] <- mean(obs.df$lnN[obs.ids])
    }
    out.df <- obs.ls[[resp]] %>% 
      mutate(yday=yday(date),
             nullGrand_lnN_lnN=mean(lnN))
  }
  return(list(yday.df=temp.df,
              obs.df=left_join(out.df, temp.df) %>% select(-yday)))
}













#' Compute F-1 and J-index across probability thresholds
#'
#' @param L.df 
#' @param prSteps 
#'
#' @return
#' @export
#'
#' @examples
compute_thresholds <- function(L.df, prSteps=0.01) {
  library(tidyverse); library(yardstick)
  pred.df <- map_dfr(seq(0, 1, by=prSteps), 
                     ~L.df %>% mutate(thresh=.x)) %>%
    mutate(pred=if_else(prA1 < thresh, "A0", "A1") %>% factor(levels=c("A0", "A1")))
  J.df <- pred.df %>%
    group_by(y, model, PCA, covSet, thresh) %>%
    j_index(truth=alert, estimate=pred, event_level="second") %>%
    rename(J=.estimate) %>%
    select(-.metric, -.estimator)
  F1.df <- pred.df %>%
    group_by(y, model, PCA, covSet, thresh) %>%
    f_meas(truth=alert, estimate=pred, event_level="second") %>%
    rename(F1=.estimate) %>%
    select(-.metric, -.estimator) 
  return(full_join(J.df, F1.df))
}














# Miscellaneous -----------------------------------------------------------




#' Shortcut for dir(..., full.names=T)
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
dirf <- function(...) {
  dir(..., full.names=T)
}










#' Adapted from butcher::axe_env for bayesian(engine='brms')
#'
#' @param x 
#' @param verbose 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
axe_env_bayesian <- function(x, verbose = FALSE, ...) {
  old <- x
  for(i in seq_along(x$fit$actions$model$spec$args)) {
    attr(x$fit$actions$model$spec$args[[i]], ".Environment") <- rlang::base_env()
  }
  for(i in seq_along(x$fit$fit$spec$args)) {
    attr(x$fit$fit$spec$args[[i]], ".Environment") <- rlang::base_env()
  }
  for(i in seq_along(x$fit$fit$spec$method$fit$args)) {
    attr(x$fit$fit$spec$method$fit$args[[i]], ".Environment") <- rlang::base_env()
  }
  attr(x$fit$actions$model$formula, ".Environment") <- rlang::base_env()

  return(x)
}




# deprecated: to delete ---------------------------------------------------



#' Extract CMEMS data to site point locations
#'
#' @param site.df 
#' @param vars 
#' @param cmems.df 
#'
#' @return
#' @export
#'
#' @examples
extract_cmems_pts <- function(site.df, vars, cmems.df) {
  library(tidyverse); library(zoo)
  
  cmems.site <- cmems.df %>% 
    filter(cmems_id %in% site_hab.df$cmems_id) %>%
    group_by(cmems_id) %>%
    mutate(across(all_of(cmems_vars), 
                  ~rollmean(.x, k=7, na.pad=T),
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
  return(cmems.site)
}




#' Extract WRF data to site point locations
#'
#' @param site.df 
#' @param site.versions 
#' @param wrf_i 
#' @param wrf.df 
#'
#' @return
#' @export
#'
#' @examples
extract_wrf_pts <- function(site.df, site.versions, wrf_i, wrf.df) {
  library(tidyverse); library(zoo)
  
  wrf.site <- wrf.df %>%
    group_by(version) %>%
    group_split() %>%
    map2_dfr(., site.versions, 
             ~.x %>% filter(wrf_id %in% site.df[[.y]])) %>%
    arrange(wrf_id, date) %>%
    group_by(wrf_id) %>%
    mutate(across(any_of(wrf_i$var), 
                  ~rollmean(.x, k=7, na.pad=T),
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
  return(wrf.site)
}




#' Extract CMEMS data to site buffers
#'
#' @param site.buffer 
#' @param vars 
#' @param cmems.df 
#'
#' @return
#' @export
#'
#' @examples
extract_cmems_buffers <- function(site.buffer, vars, cmems.df) {
  
  library(tidyverse); library(zoo)
  
  cmems.buffer <- expand_grid(siteid=unique(site.buffer$siteid),
                              quadrant=unique(site.buffer$quadrant),
                              date=unique(cmems.df$date)) %>%
    bind_cols(as_tibble(setNames(map(cmems_vars, ~NA_real_), cmems_vars)))
  
  cmems_id.ls <- map(site.buffer$cmems_id, ~.x)
  cmems.df$date_id <- match(cmems.df$date, unique(cmems.df$date)) 
  cmems_dates.ls <- map(1:n_distinct(cmems.df$date_id), ~which(cmems.df$date_id==.x))
  
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
                  ~rollmean(.x, k=7, na.pad=T),
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
    mutate(across(where(is.numeric), na.aggregate)) %>%
    ungroup
  
  return(cmems.buffer)
}





#' Extract WRF data to site buffers
#'
#' @param site.buffer 
#' @param wrf_i 
#' @param wrf.df 
#'
#' @return
#' @export
#'
#' @examples
extract_wrf_buffers <- function(site.buffer, wrf_i, wrf.df) {
  library(tidyverse); library(sf); library(zoo)
  
  wrf.buffer <- expand_grid(siteid=unique(site.buffer$siteid),
                            quadrant=unique(site.buffer$quadrant),
                            date=unique(wrf.df$date)) %>%
    mutate(version=1 + (date >= "2019-04-01")) %>%
    bind_cols(as_tibble(setNames(map(wrf_i$vars, ~NA_real_), wrf_i$vars)))
  
  wrf_id.ls <- list(v1=map(site.buffer$wrf_id.1, ~.x),
                    v2=map(site.buffer$wrf_id.2, ~.x))
  wrf.df$date_id <- match(wrf.df$date, unique(wrf.df$date)) 
  wrf_dates.ls <- map(1:n_distinct(wrf.df$date_id), ~which(wrf.df$date_id==.x))
  
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
                  ~rollmean(.x, k=7, na.pad=T),
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
    mutate(across(where(is.numeric), na.aggregate)) %>%
    ungroup
  return(wrf.buffer)
}




