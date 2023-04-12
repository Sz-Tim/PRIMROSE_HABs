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
  
  for(i in seq_along(wrf_dates)) {
    wrf_i.i <- filter(wrf_i, date_0==wrf_dates[i])
    nc_f.i <- map(wrf_i.i$fname, ~glue("{wrf.dir}/{.x}")) %>% 
      set_names(wrf_i.i$res)
    d_i <- nest_WRF_domains(nc_f.i)
    nc.ls <- map(nc_f.i, nc_open)
    time.ls <- map(nc.ls, 
                   ~tibble(Times=ncvar_get(.x, "Times")) %>%
                     mutate(Time.dt=as_datetime(str_replace(Times, "_", " "))))
    var.df <- map2(nc.ls, time.ls,
                   ~expand_grid(time=.y$Time.dt,
                                lat_i=1:(.x$dim$south_north$len),
                                lon_i=1:(.x$dim$west_east$len)) %>%
                     mutate(date=date(time),
                            U=c(ncvar_get(.x, "U10")),
                            V=c(ncvar_get(.x, "V10")),
                            Shortwave=c(ncvar_get(.x, "Shortwave")),
                            Precipitation=c(ncvar_get(.x, "Precipitation")),
                            sst=c(ncvar_get(.x, "sst")))) %>%
      map2_dfr(., d_i,
               ~.x %>% 
                 group_by(time) %>%
                 mutate(i=row_number()) %>%
                 ungroup %>%
                 right_join(., .y %>% select(i, elev, res, lon, lat)) %>%
                 group_by(i, res, lat, lon, date) %>%
                 summarise(UV_90=q90(sqrt(U^2 + V^2)),
                           UV_mn=mean(sqrt(U^2 + V^2)),
                           UV_mnDir=atan2(mean(V), mean(U)),
                           Shortwave=sum(Shortwave),
                           Precip=sum(Precipitation),
                           sst=mean(sst)) %>%
                 ungroup)
    saveRDS(var.df, glue("{out.dir}/wrf/wrf_{str_pad(i, 4, 'left', '0')}.rds"))
  }
  
  
  
  
  
  out.ls <- vector("list", n_distinct(sampling.df$wrf_i))
  for(i in unique(sampling.df$wrf_i)) {
    samp_i <- sampling.df %>% filter(wrf_i==i) %>%
      mutate(wrf_coord=st_nearest_feature(., coord.wrf),
             wrf_row=coord.wrf$row[wrf_coord],
             wrf_col=coord.wrf$col[wrf_coord])
    
  
    
    for(j in 1:nrow(samp_i)) {
      dims_ij <- list(lon=samp_i$wrf_row[j],
                      lat=samp_i$wrf_col[j],
                      time=which(times.df$date == samp_i$date[j]))
      U_ij <- U[dims_ij$lon, dims_ij$lat, dims_ij$time]
      V_ij <- V[dims_ij$lon, dims_ij$lat, dims_ij$time]
      var_ij[[j]] <- tibble(
        obs.id=samp_i$obs.id[j],
        site.id=samp_i$site.id[j],
        date=samp_i$date[j],
        U_mn=mean(U_ij),
        V_mn=mean(V_ij),
        UV_speed=q90(sqrt(U_ij^2 + V_ij^2)),
        UV_mnDir=mean(atan2(V_ij, U_ij)),
        Shortwave=sum(Shortwave[dims_ij$lon, dims_ij$lat, dims_ij$time]),
        Precip=sum(Precip[dims_ij$lon, dims_ij$lat, dims_ij$time]),
        sst=mean(sst[dims_ij$lon, dims_ij$lat, dims_ij$time])
      )
    }
    out.ls[[i]] <- do.call('rbind', var_ij)
  }
}






nest_WRF_domains <- function(path.ls) {
  library(tidyverse); library(ncdf4); library(sf)
  nDomain <- length(path.ls)
  nc.ls <- map(path.ls, nc_open)
  lon <- map(nc.ls, ~ncvar_get(.x, "XLONG"))
  lat <- map(nc.ls, ~ncvar_get(.x, "XLAT"))
  elev <- map(nc.ls, ~ncvar_get(.x, "HGT"))
  walk(nc.ls, nc_close)
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





get_shortestPaths <- function(ocean.path, site.df) {
  library(tidyverse); library(sf); library(glue);
  library(raster); library(gdistance)
  
  # adapted from:
  # https://agrdatasci.github.io/gdistance/reference/index.html
  
  # load ocean raster and calculate transition matrix
  mesh.r <- raster(ocean.path)
  mesh.tmx <- transition(mesh.r, mean, 16)
  mesh.tmx <- geoCorrection(mesh.tmx)
  
  # locate sites in mesh
  set_ll_warn(TRUE)
  site.spdf <- SpatialPointsDataFrame(site.df[,c("lon", "lat")],
                                      data=site.df[,"siteid"], 
                                      proj4string=CRS("+init=epsg:27700")) %>%
    points2nearestcell(., mesh.r) %>%
    as.data.frame
  
  # find pairwise shortest paths within ocean
  paths.ls <- list()
  for(i in 1:nrow(site.spdf)) {
    paths.ls[[i]] <- shortestPath(mesh.tmx, 
                                  as.matrix(site.spdf[i,2:3]),
                                  as.matrix(site.spdf[-i,2:3]),
                                  output="SpatialLines") %>%
      st_as_sf %>%
      mutate(origin=site.spdf$site.id[i],
             destination=site.spdf$site.id[-i],
             distance=st_length(.))
  }
  
  dist.df <- map_dfr(1:nrow(site.spdf),
                     ~shortestPath(mesh.tmx,
                                   as.matrix(site.spdf[.x,2:3]),
                                   as.matrix(site.spdf[-.x,2:3]),
                                   output="SpatialLines") %>%
                       st_as_sf %>%
                       mutate(origin=site.spdf$site.id[.x],
                              destination=site.spdf$site.id[-.x],
                              distance=st_length(.))) %>%
    st_drop_geometry()
  
  return(dist.df)
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
