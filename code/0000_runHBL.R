# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Model periodic re-fits



# setup -------------------------------------------------------------------
pkgs <- c("tidyverse", "glue", "tidymodels", 
          "nnet", "randomForest", "glmnet", "xgboost", "earth",
          "brms", "bayesian", "doParallel", "foreach", "butcher")
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
  arrange(y, id) 
cores_per_model <- 3
n_spp_parallel <- 18

# i <- 14
registerDoParallel(n_spp_parallel)
foreach(i=1:nrow(covSet.df)) %dopar% {
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  base.dir <- "out/0_init" 
  
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
      paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X"),
      paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X")
    ),
    hab=c(outer(filter(y_i, type=="hab")$abbr, c("lnNAvg", "prA"), "paste0"))
  )
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
  saveRDS(obs.split, glue("data/0_init/{y}_{covSet}_dataSplit.rds"))
  obs.train <- training(obs.split)
  obs.test <- testing(obs.split)
  
  cat("Starting", covSet, "for", y, ":", as.character(Sys.time()), "\n",
      file=glue("out/logs/d{d}_{y}_HBL.log"))
  
  
  
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
    ~list(HBL=make_HB_formula(.x, c(covs$main, covs$interact)),
          HBL_PCA=make_HB_formula(.x, covsPCA),
          HB_vars=formula(glue("{.x} ~ yday + lon + lat + siteid +",
                               "{paste(unlist(covs), collapse='+')}")),
          HB_vars_PCA=formula(glue("{.x} ~ yday + lon + lat + siteid +",
                                   "{paste(covsPCA, collapse='+')}"))
    )
  )
  
  # priors
  priStr <- switch(1,
                   "1"=list(r1=0.2, r2=2, hs1=0.5, hs2=0.6, b=0.75, de=0.3, i=1),
                   "2"=list(r1=0.1, r2=2, hs1=3, hs2=0.2, b=0.2, de=0.1, i=2),
                   "3"=list(r1=0.1, r2=1, hs1=5, hs2=0.5, b=0.5, de=0.05, i=3)
  )
  priors <- map(responses, ~list(HBL=make_HB_priors(priStr, "HBL", .x, covs),
                                 HBL_PCA=make_HB_priors(priStr, "HBL", .x, covsPCA, PCA=T)))
  
  # tuning controls
  HB.i <- list(
    iter=1500,
    warmup=1000,
    chains=cores_per_model,
    cores=cores_per_model,
    ctrl=list(adapt_delta=0.8, max_treedepth=10),
    prior_i=priStr$i
  )
  
  
  
  # . train: fit models -----------------------------------------------------
  
  for(r in responses) {
    set.seed(1003)
    folds <- vfold_cv(d.y$train[[r]], strata=r)
    set.seed(1003)
    foldsPCA <- vfold_cv(dPCA.y$train[[r]], strata=r)
    fit_model("HBL", r, form.ls, dPCA.y$train, HB.i, priors, fit.dir, y, "_PCA")
    fit_model("HBL", r, form.ls, d.y$train, HB.i, priors, fit.dir, y)
  }
  
  
  
  # . cross validation by year ----------------------------------------------
  
  # Only necessary for Bayesian models; ML models use CV in fitting process
  for(r in responses) {
    set.seed(1003)
    folds <- vfold_cv(d.y$train[[r]], strata=r)
    set.seed(1003)
    foldsPCA <- vfold_cv(dPCA.y$train[[r]], strata=r)
    run_Bayes_CV("HBL", foldsPCA, cv.dir, y, y_i.i, r, form.ls, HB.i, priors, PCA=T)
    run_Bayes_CV("HBL", folds, cv.dir, y, y_i.i, r, form.ls, HB.i, priors)
  }
  cat("Finished", covSet, "for", y, ":", as.character(Sys.time()), "\n")
}


closeAllConnections()




