# Calculation of Accumulated Local Effects

# NB: The structure is strange because I'm calculating HBL on my desktop locally
# because they were fit using brms 2.19.0 and R 4.3.0, and do *not* work with
# newer versions. R was updated to 4.4 on Salmon, and given time constraints 
# with the revisions, this was the best solution. Toggle the 'mods' call to
# switch back and forth and everything should work if run on the correct machine.
# Future revisions will use the renv package to manage these issues.

# setup -------------------------------------------------------------------

mods <- c("Ridge", "MARS", "RF", "NN", "Boost")
mods <- "HBL1"

pkgs <- c("tidyverse", "tidymodels", "brms",
          "bayesian", "doParallel", "foreach", "butcher", "ingredients")
if(any(mods != "HBL1")) {
  pkgs <- c(pkgs, "nnet", "randomForest", "glmnet", "xgboost", "earth")
}
pkg_dir <- "~/R/old_versions"
lapply(pkgs, library, character.only=T, lib.loc=pkg_dir)
library(glue)
source("code/00_fn.R")

y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv", show_col_types=F) |> 
                   arrange(abbr) |> mutate(type="tox")) |>
  filter(! abbr %in% c("AZP", "YTX"))
y_resp <- filter(y_i, abbr!="Prli")$abbr

covSet.df <- expand_grid(y=y_resp,
                         Avg=c(0,1), 
                         Xf=c(0,1),
                         XN=c(0,1),
                         Del=c(0,1)) |>
  group_by(y) |>
  mutate(id=row_number(),
         f=glue("{id}-Avg{Avg}_Xf{Xf}_XN{XN}_Del{Del}")) |>
  ungroup()

if(all(mods=="HBL1")) {
  n_spp_parallel <- 2
  base.dir <- "D:/PRIMROSE_HABs_submitted/PRIMROSE_HABs/out/0_init"
} else {
  n_spp_parallel <- 4
  base.dir <- "out/0_init" 
}



accumulated_dependence_revised <- function(x,
                                           data,
                                           predict_function = predict,
                                           label = class(x)[1],
                                           variables = NULL,
                                           N = 500,
                                           variable_splits = NULL,
                                           grid_points = 101,
                                           ...,
                                           variable_type = "numerical") {
  if (!is.null(N) && N < nrow(data)) {
    # sample N points
    ndata <- data[sample(1:nrow(data), N), , drop = FALSE]
  } else {
    ndata <- data
  }
  
  cp <- ceteris_paribus(x,
                        data,
                        predict_function = predict_function,
                        new_observation = ndata,
                        variables = variables,
                        grid_points = grid_points,
                        variable_splits = variable_splits,
                        label = label, ...)
  
  ALE <- aggregate_profiles(cp, variables = variables, type = "accumulated", variable_type = variable_type, ...)
  return(list(CP=cp, ALE=ALE))
}




# Model predictions -------------------------------------------------------

registerDoParallel(n_spp_parallel)
foreach(i=(1:nrow(covSet.df))) %dopar% {
  
  lapply(pkgs, library, character.only=T, lib.loc=pkg_dir)
  library(glue)
  source("code/00_fn.R")
  if(all(mods=="HBL1")) {
    base.dir <- "D:/PRIMROSE_HABs_submitted/PRIMROSE_HABs/out/0_init"
  } else {
    base.dir <- "out/0_init" 
  }
  
  covSet <- covSet.df$f[i]
  d <- covSet.df$id[i]
  
  fit.dir <- glue("{base.dir}/model_fits/{covSet}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("{base.dir}/ensembles/")
  out.dir <- glue("{base.dir}/compiled/{covSet}/")
  dir.create(ens.dir, recursive=T, showWarnings=F)
  dir.create(cv.dir, recursive=T, showWarnings=F)
  dir.create(out.dir, recursive=T, showWarnings=F)
  dir.create(glue("out/0_init/compiled/{covSet}/"))
  
  y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) |> 
                     arrange(abbr) |> mutate(type="hab"),
                   read_csv("data/i_tox.csv", show_col_types=F) |> 
                     arrange(abbr) |> mutate(type="tox")) |>
    filter(! abbr %in% c("AZP", "YTX"))
  
  y <- covSet.df$y[i]
  y_i.i <- filter(y_i, abbr==y)
  cat("Starting", y, d, "at", as.character(Sys.time()), "\n", 
      file=glue("{base.dir}/logs/{y}_{str_pad(d, 2, 'left', '0')}.txt"))
  
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
      "lnNPrevYr", "lnNAvgPrevYr", "prAlertPrevYr", "prAlertAvgPrevYr",
      col_cmems, col_wrf
    ),
    interact=c(
      paste("UWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
      paste("VWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
      paste("UWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
      paste("VWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
      paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X"),
      paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[NS]", col_wrf, value=T), sep="X"),
      paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[NS]", col_wrf, value=T), sep="X"),
      paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X")
    ),
    hab=c(outer(filter(y_i, type=="hab")$abbr, c("lnNAvg", "prA"), "paste0"))
  )
  all_covs$interact <- paste("lnNWt1", c(all_covs$main[-2]), sep="X")
  
  covs_exclude <- get_excluded_cov_regex(covSet)
  
  obs.ls <- map_dfr(dirf("data/0_init", "data_.*_all.rds"), readRDS) |>
    filter(y %in% y_i$abbr) |>
    filter(y != "Prli") |>
    filter(year(date) < 2023) |>
    select(all_of(col_metadata), all_of(col_resp),
           "alert1", "alert2", any_of(unname(unlist(all_covs)))) |>
    mutate(across(starts_with("alert"), ~factor(.x)),
           across(starts_with("tl"), ~factor(.x, ordered=T))) |>
    group_by(y, obsid) |>
    slice_head(n=1) |>
    filter(y==y_i.i$abbr) |>
    select(where(~any(!is.na(.x)))) |>
    na.omit()
  
  set.seed(1003)
  obs.split <- group_initial_split(obs.ls, group=year)
  obs.train <- training(obs.split)
  obs.test <- testing(obs.split)
 
  
  
  
  
  # . prep ------------------------------------------------------------------
  
  responses <- c(alert="alert", tl="tl")[1]
  prep.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude))
  d.y <- list(train=map(prep.ls, ~bake(.x, obs.train)),
              test=map(prep.ls, ~bake(.x, obs.test)))
  covs <- filter_corr_covs(all_covs, d.y, covs_exclude)
  covs_ALE <- c(covs$main, "yday", "lat", "lon") |>
    grep("alert", x=_, value=T, invert=T)
  if(y_i.i$type == "tox") {
    covs_ALE <- c(covs_ALE, covs$hab)
  }
  
  for(p in 1:2) {
    PCA <- c(TRUE, FALSE)[p]
    if(PCA) {
      prep_p <- prep_recipe(obs.train, "alert", covs_exclude, TRUE)
    } else {
      prep_p <- prep_recipe(obs.train, "alert", covs_exclude)
    }
    
    for(m in seq_along(mods)) {
      mod_m <- mods[m]
      if(file.exists(glue("out/0_init/compiled/{covSet}/{y}_adp_{mod_m}{ifelse(PCA, '_PCA', '')}.rds"))) next
      
      get_predictions_i <- function(fit, data_orig, mod=mod_m, resp="alert", y_i.i=y, 
                                    PCA=PCA, prep_p_=prep_p) {
        library(tidyverse, lib.loc="~/R/old_versions"); 
        library(glue); 
        library(tidymodels, lib.loc="~/R/old_versions"); 
        library(brms, lib.loc="~/R/old_versions")
        
        d.df <- bake(prep_p_, data_orig)
        
        if(grepl("HB", mod)) {
          pred <- parsnip::extract_fit_engine(fit) |>
            posterior_epred(d.df, allow_new_levels=T) |>
            summarise_post_preds(resp, y_i.i)
        } else {
          pred_type <- ifelse(resp=="lnN", "raw", "prob")
          preds <- predict(fit, d.df, pred_type)
          pred <- summarise_ML_preds(preds, resp, y_i.i)
        }
        return(c(pred))
      }
      
      set.seed(1003)
      mod_path <- glue("{base.dir}/model_fits/{covSet}/{y}_alert_{mod_m}{ifelse(PCA, '_PCA', '')}.rds")
      if(grepl("HB", mod_m)) {
        adp <- accumulated_dependence_revised(readRDS(mod_path), 
                                              as.data.frame(obs.ls), 
                                              covs_ALE, 
                                              label=mod_m, 
                                              predict_function=get_predictions_i)
      } else {
        adp <- accumulated_dependence_revised(readRDS(mod_path), 
                                              as.data.frame(obs.ls), 
                                              covs_ALE, 
                                              label=mod_m, 
                                              predict_function=get_predictions_i)
      }
      adp_df <- adp$ALE |>
        set_names(c("variable", "model", "x", "yhat", "id")) |>
        select(variable, x, yhat)
      names(adp_df)[3] <- paste0("d", d, ".", ifelse(PCA, "PCA.", ""), mod_m, "_alert_A1")
      saveRDS(adp_df, glue("out/0_init/compiled/{covSet}/{y}_adp_{mod_m}{ifelse(PCA, '_PCA', '')}.rds"))
      rm(adp_df); gc()
      cp_df <- adp$CP |>
        group_by(`_vname_`, `_ids_`) |>
        mutate(cp_x_id=row_number()) |>
        ungroup()
      if(mod_m=="Ridge" & PCA) {
        saveRDS(cp_df |> select(-`_yhat_`, -`_label_`), 
                glue("out/0_init/compiled/{covSet}/{y}_cp_data.rds"))
        saveRDS(attr(adp$CP, "observations"), 
                glue("out/0_init/compiled/{covSet}/{y}_cp_data_attr.rds"))
      }
      cp_df <- cp_df |>
        select(`_vname_`, `_ids_`, cp_x_id, `_yhat_`)
      names(cp_df)[4] <- paste0("d", d, ".", ifelse(PCA, "PCA.", ""), mod_m, "_alert_A1")
      saveRDS(cp_df, glue("out/0_init/compiled/{covSet}/{y}_cp_{mod_m}{ifelse(PCA, '_PCA', '')}.rds"))
      rm(cp_df); gc()
      cat("  Finished", m, "PCA:", PCA, "at", as.character(Sys.time()), "\n", 
          file=glue("{base.dir}/logs/{y}_{str_pad(d, 2, 'left', '0')}.txt"), append=T)
    }
  }
  cat("Finished", y, d, "at", as.character(Sys.time()), "\n", 
      file=glue("{base.dir}/logs/{y}_{str_pad(d, 2, 'left', '0')}.txt"), append=T)
}

closeAllConnections()
