
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
  # download nc files
  command <- paste("copernicusmarine subset -i", " subset -i", i.df$ID_toolbox[i],
                   "-x", bbox$xmin, "-X", bbox$xmax,
                   "-y", bbox$ymin, "-Y", bbox$ymax,
                   "-t", dateRng[1]-nDays_buffer, "-T", dateRng[2]+nDays_buffer,
                   " --variable", i.df$var[i],
                   "-o temp", "-f", i.df$fname[i],
                   "--force-download --overwrite-metadata-cache")
  system(command, intern=TRUE)
}

# Process .nc files
for(i in 1:nrow(i.df)) {
  # load .nc file
  nc <- nc_open(paste0(out.dir, i.df$fname[i]))
  
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

# Remove temporary .nc files
file.remove(i.df$fname)
