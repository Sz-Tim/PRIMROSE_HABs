# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Functions





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
    if(i.df$var[i]=="spco2") {
      var.start <- c(lon_i[1], lat_i[1], time_i[1])
      var.count <- c(length(lon_i), length(lat_i), length(time_i)) 
    } else {
      var.start <- c(lon_i[1], lat_i[1], 1, time_i[1])
      var.count <- c(length(lon_i), length(lat_i), 1, length(time_i)) 
    }
    
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








get_WRF <- function(wrf.dir, nDays_buffer, dateRng, out.dir) {
  library(tidyverse); library(ncdf4); library(lubridate); library(glue)
  dir.create(glue("{out.dir}/wrf"), showWarnings=F)
  
  # metadata for all WRF files within timespan
  wrf_i <- tibble(fname=dir(wrf.dir, ".nc$", recursive=T)) %>%
    mutate(res=str_sub(fname, -6, -4),
           day_0=str_sub(fname, 23, 24),
           month_0=str_sub(fname, 21, 22),
           year_0=str_sub(fname, 17, 20),
           day_1=str_sub(fname, 28, 29),
           month_1=str_sub(fname, 26, 27),
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
    nc_f.i <- map(wrf_i.i$fname, ~glue("{wrf.dir}/{.x}")) %>% 
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
      j.fname <- glue("{out.dir}/wrf/wrf_{var.ls[[j]]$date[1]}_{var.ls[[j]]$res[1]}.rds")
      if(!file.exists(j.fname)) {
        var.ls[[j]] %>% 
          group_by(i, lat_i, lon_i, date) %>%
          summarise(UV_90=q90(sqrt(U^2 + V^2)),
                    UV_mn=mean(sqrt(U^2 + V^2)),
                    UV_mnDir=atan2(mean(V), mean(U)),
                    Shortwave=sum(Shortwave),
                    Precip=sum(Precipitation),
                    sst=mean(sst)) %>%
          ungroup %>%
          saveRDS(j.fname)
      }
    }
    # save domain extents IF all three domains are represented
    if(length(nc_f.i)==3 & all(c("d01", "d02", "d03") %in% names(nc_f.i))) {
      nest_WRF_domains(nc.ls) %>%
        walk(~saveRDS(.x, glue("{out.dir}/wrf/wrfDomains_{wrf_dates[i]}_{.x$res[1]}.rds")))
      
    }
    walk(nc.ls, nc_close)
  }
}






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





get_shortestPaths <- function(ocean.path, site.df, transMx.path=NULL, recalc_transMx=T) {
  library(tidyverse); library(sf); library(glue);
  library(raster); library(gdistance)
  
  # adapted from:
  # https://agrdatasci.github.io/gdistance/reference/index.html
  
  # load ocean raster and calculate transition matrix
  mesh.r <- raster(ocean.path)
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
  
  # find pairwise shortest paths within ocean
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





get_fetch <- function(site.df, fetch.path, small=T, buffer=5e2, fun=mean) {
  library(tidyverse); library(sf)
  site.df %>%
    st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F) %>%
    mutate(fetch=raster::extract(raster(fetch.path), ., 
                                 small=small, buffer=buffer, fun=fun)) %>%
    st_drop_geometry()
}





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





get_trafficLights <- function(hab.df, N, tl.df) {
  library(tidyverse)
  hab.df %>%
    rowwise() %>%
    mutate(tl=tl.df$tl[max(which(sp==tl.df$sp & {{N}} >= tl.df$min_ge))],
           cat=tl.df$N_ord[max(which(sp==tl.df$sp & {{N}} >= tl.df$min_ge))],
           alert=tl.df$alert[max(which(sp==tl.df$sp & {{N}} >= tl.df$min_ge))]) %>%
    ungroup %>%
    mutate(alert=c("0_none", "1_warn", "2_alert")[alert+1])
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





# modified from astsa::trend
detrend_loess <- function (x, y, span=0.75, robust=TRUE) {
  if(length(y) < 10) {
    return(y)
  }
  fam = ifelse(robust, "symmetric", "gaussian")
  lo = stats::predict(stats::loess(y ~ x, span=span, family=fam), se = F)
  return(c(y - lo))
}
