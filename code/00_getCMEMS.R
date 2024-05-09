
# Everything updated to work with new Copernicus data access methods
# This script needs to be run inside a virtual environment to be able to 
# load the (mandatory) copernicus marine toolbox. 
# Consequently, get_CMEMS() creates a temporary .RData file with the necessary
# objects and then calls a bash script which starts the virtual environment 
# and runs this script. This is all much more complicated since it was really
# built under the assumption that everyone would use Python....

library(tidyverse); library(terra); library(lubridate); library(glue)

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

for(i in unique(i.df$var)) {
  nc_i <- dir("temp", paste0(i, ".*nc"), full.names=T)
  if(length(nc_i) == 1) {
    nc <- rast(nc_i)
    
  } else {
    nc_F <- rast(nc_i |> grep("Forecast", x=_, value=T)) 
    nc_R <- rast(nc_i |> grep("Reanalysis", x=_, value=T)) |>
      project(nc_F)
    nc <- mergeTime(sds(nc_R, nc_F), mean) 
  }
  nc_LU <- crds(nc) |> as_tibble() |>
    rename(lon=x, lat=y) |>
    mutate(cmems_id=row_number())
  saveRDS(nc_LU, glue("{out.dir}/coords_{i}.rds"))
  nc_df <- nc_LU |> select(cmems_id) |>
    bind_cols(values(nc, dataframe=T, na.rm=T) |>
                setNames(time(nc))) |>
    pivot_longer(-1, names_to="date", values_to=i)
  saveRDS(nc_df, glue("{out.dir}/cmems_{i}.rds"))
}
