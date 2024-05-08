
# Everything updated to work with new Copernicus data access methods
# This script needs to be run inside a virtual environment to be able to 
# load the (mandatory) copernicus marine toolbox. 
# Consequently, get_CMEMS() creates a temporary .RData file with the necessary
# objects and then calls a bash script which starts the virtual environment 
# and runs this script. This is all much more complicated since it was really
# built under the assumption that everyone would use Python....

library(tidyverse); library(ncdf4); library(lubridate); library(glue)

load("temp/get_CMEMS.RData")
i.df <- i.df |>
  mutate(fname=glue("cmems_{var}_{source}.nc"))
for(i in 1:nrow(i.df)) {
  if(grepl("Forecast", i.df$source[i])) {
    dates <- c(pmax(ymd("2021-09-01"), dateRng[1]-nDays_buffer),
               pmin(today()+3, dateRng[2]+nDays_buffer))
  } else {
    dates <- c(pmax(ymd("1993-01-01"), dateRng[1]-nDays_buffer),
               pmin(ymd("2024-01-01"), dateRng[2]+nDays_buffer))
  }
  # download nc files
  command <- paste("copernicusmarine subset -i", i.df$ID_toolbox[i],
                   "-x", bbox$xmin, "-X", bbox$xmax,
                   "-y", bbox$ymin, "-Y", bbox$ymax,
                   "-z", 0, "-Z", 0,
                   "-t", dates[1], "-T", dates[2],
                   " -v", i.df$var[i],
                   "-o temp", "-f", i.df$fname[i],
                   "--force-download --overwrite-metadata-cache")
  system(command, intern=TRUE)
}

nc_f <- dir("temp", "cmems_.*nc")
# Process .nc files
for(i in 1:length(nc_f)) {
  # load .nc file
  nc <- nc_open(paste0("temp/", nc_f[i]))
  
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
    mutate(source=str_sub(str_split_fixed(nc_f[i], "_", 3)[,3], 1, -4),
           vals=c(nc.var)) |>
    rename_with(~gsub("vals", str_split_fixed(nc_f[i], "_", 3)[,2], .x)) |>
    saveRDS(glue("{out.dir}/{str_replace(nc_f[i], 'nc', 'rds')}"))
}

# Remove temporary .nc files
file.remove(nc_f)
