# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Model periodic re-fits


# setup -------------------------------------------------------------------
pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", 
          "nnet", "randomForest", "glmnet", "xgboost", "earth",
          "brms", "bayesian", "doParallel", "butcher")
lapply(pkgs, library, character.only=T)
source("code/00_fn.R")

y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) %>% 
                   arrange(abbr) %>% mutate(type="hab"),
                 read_csv("data/i_tox.csv", show_col_types=F) %>% 
                   arrange(abbr) %>% mutate(type="tox")) %>%
  filter(! abbr %in% c("AZP", "YTX"))
y_resp <- filter(y_i, abbr!="Prli")$abbr

covSet.df <- expand_grid(y=y_resp,
                         Avg=c(0,1), 
                         Xf=c(0,1),
                         XN=c(0,1),
                         Del=c(0,1)) %>%
  group_by(y) %>%
  mutate(id=row_number(),
         f=glue("{id}-Avg{Avg}_Xf{Xf}_XN{XN}_Del{Del}")) %>%
  ungroup %>%
  arrange(y, id) |>
  filter(Xf==1 & Avg==1)
n_spp_parallel <- 1
cores_per_model <- 12


# registerDoParallel(cores_per_model)
for(i in 1:nrow(covSet.df)) {
# registerDoParallel(n_spp_parallel)
# foreach(i=1:nrow(covSet.df)) %dopar% {
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  base.dir <- "out/0_init_redo" 
  
  covSet <- covSet.df$f[i]
  d <- covSet.df$id[i]
  fit.dir <- glue("{base.dir}/model_fits/{covSet}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("{base.dir}/ensembles/")
  out.dir <- glue("{base.dir}/compiled/{covSet}/")
  dir.create(ens.dir, recursive=T, showWarnings=F)
  dir.create(cv.dir, recursive=T, showWarnings=F)
  dir.create(out.dir, recursive=T, showWarnings=F)
  dir.create(glue("{fit.dir}/vi/"), recursive=T, showWarnings=F)
  
  y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) %>% 
                     arrange(abbr) %>% mutate(type="hab"),
                   read_csv("data/i_tox.csv", show_col_types=F) %>% 
                     arrange(abbr) %>% mutate(type="tox")) %>%
    filter(! abbr %in% c("AZP", "YTX"))
  
  y <- covSet.df$y[i]
  y_i.i <- filter(y_i, abbr==y)
  
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
  # all_covs$interact <- paste("lnNWt1", c(all_covs$main[-2]), sep="X")
  all_covs$interact <- c(all_covs$interact,
                         paste("lnNWt1", c(all_covs$main[-2]), sep="X"))
  
  covs_exclude <- get_excluded_cov_regex(covSet)
  
  obs.ls <- map_dfr(dirf("data/0_init", "data_.*_all.rds"), readRDS) %>%
    filter(y %in% y_i$abbr) %>%
    filter(y != "Prli") %>%
    filter(year(date) < 2023) %>%
    select(all_of(col_metadata), all_of(col_resp),
           "alert1", "alert2", any_of(unname(unlist(all_covs)))) %>%
    mutate(across(starts_with("alert"), ~factor(.x)),
           across(starts_with("tl"), ~factor(.x, ordered=T))) %>%
    group_by(y, obsid) %>%
    slice_head(n=1) %>%
    filter(y==y_i.i$abbr) %>%
    select(where(~any(!is.na(.x)))) %>%
    na.omit()
  
  set.seed(1003)
  obs.split <- group_initial_split(obs.ls, group=year)
  saveRDS(obs.split, glue("data/0_init_redo/{y}_{covSet}_dataSplit.rds"))
  obs.train <- training(obs.split)
  obs.test <- testing(obs.split)
  
  cat("Starting", covSet, "for", y, ":", as.character(Sys.time()), "\n",
      file=glue("out/logs/d{d}_{y}_ML.log"))
  
  
  
  # . prep ------------------------------------------------------------------
  
  responses <- c(alert="alert", tl="tl")[1]
  prep.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude))
  d.y <- list(train=map(prep.ls, ~bake(.x, obs.train)),
              test=map(prep.ls, ~bake(.x, obs.test)))
  prepPCA.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude, TRUE))
  dPCA.y <- list(train=map(prepPCA.ls, ~bake(.x, obs.train)),
                 test=map(prepPCA.ls, ~bake(.x, obs.test)))
  d.split <- dPCA.split <- map(responses, ~obs.split)
  for(r in responses) {
    d.split[[r]]$data <- d.split[[r]]$data %>% select(y, obsid) %>%
      left_join(bind_rows(d.y$train[[r]], d.y$test[[r]]), by=c("y", "obsid"))
    dPCA.split[[r]]$data <- dPCA.split[[r]]$data %>% select(y, obsid) %>%
      left_join(bind_rows(dPCA.y$train[[r]], dPCA.y$test[[r]]), by=c("y", "obsid"))
  }
  
  covs <- filter_corr_covs(all_covs, d.y, covs_exclude)
  covsPCA <- names(dPCA.y$train[[1]] %>% select(starts_with("PC")))
  
  # formulas
  form.ls <- map(
    responses,
    ~list(ML=formula(glue("{.x} ~ {paste(unlist(covs), collapse='+')}")),
          ML_PCA=formula(glue("{.x} ~ {paste(covsPCA, collapse='+')}"))
    )
  )
  
  # tuning controls
  n_tuneVal <- list(Ridge=1e3,
                    MARS=1e3,
                    RF=1e2, 
                    NN=1e2,
                    Boost=1e2)
  
  
  
  # . train: fit models -----------------------------------------------------
  
  cl <- makeCluster(cores_per_model)
  registerDoParallel(cl)
  for(r in responses) {
    set.seed(1003)
    folds <- vfold_cv(d.y$train[[r]], strata=r)
    set.seed(1003)
    foldsPCA <- vfold_cv(dPCA.y$train[[r]], strata=r)
    # fit_model("Ridge", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    # fit_model("Ridge", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    # fit_model("MARS", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    # fit_model("MARS", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("RF", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    fit_model("RF", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    # fit_model("NN", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    # fit_model("NN", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    # fit_model("Boost", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    # fit_model("Boost", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
  }
  stopCluster(cl)
  closeAllConnections()
  
  cat("Finished", covSet, "for", y, ":", as.character(Sys.time()), "\n")
}

closeAllConnections()















