# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model fits



# setup -------------------------------------------------------------------
pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", "nnet", "randomForest", "glmnet", 
          "brms", "bayesian", "doParallel", "foreach", "butcher")
lapply(pkgs, library, character.only=T)
source("code/00_fn.R")

covSet <- c("1-noDtDeltaX", "2-noX", "3-noDtDelta", "4-noDt", "5-full")[1]
cores_per_model <- 5
n_spp_parallel <- 5
fit.dir <- glue("out/model_fits/{covSet}/")
cv.dir <- glue("{fit.dir}/cv/")
ens.dir <- "out/model_fits/ensembles/"
out.dir <- glue("out/compiled/{covSet}/")
dir.create(fit.dir, recursive=T, showWarnings=F)
dir.create(cv.dir, recursive=T, showWarnings=F)
dir.create(ens.dir, recursive=T, showWarnings=F)
dir.create(out.dir, recursive=T, showWarnings=F)

y_i <- bind_rows(read_csv("data/i_hab.csv") %>% arrange(abbr) %>% mutate(type="hab"),
                 read_csv("data/i_tox.csv") %>% arrange(abbr) %>% mutate(type="tox")) %>%
  filter(! abbr %in% c("ASP", "AZP", "YTX"))

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

covs_exclude <- get_excluded_cov_regex(covSet)

obs.ls <- map_dfr(dirf("data/0_init", "data_.*_all.rds"), readRDS) %>%
  filter(y %in% y_i$abbr) %>%
  filter(year(date) < 2023) %>%
  select(all_of(col_metadata), all_of(col_resp), 
         "alert1", "alert2", any_of(unname(unlist(all_covs)))) %>%
  mutate(across(starts_with("alert"), ~factor(.x)),
         across(starts_with("tl"), ~factor(.x, ordered=T))) %>%
  group_by(y, obsid) %>%
  slice_head(n=1) %>%
  group_by(y) %>%
  group_split() %>%
  map(~.x %>% select(where(~any(!is.na(.x)))) %>% na.omit)



# runs by taxon -----------------------------------------------------------

registerDoParallel(n_spp_parallel)
foreach(i=seq_along(obs.ls),
        .export=c("all_covs", "obs.ls", "covSet", "covs_exclude", 
                  "fit.dir", "out.dir", "cv.dir", "ens.dir", "y_i")
) %dopar% {
  
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  y <- obs.ls[[i]]$y[1]
  y_i.i <- filter(y_i, abbr==y)
  set.seed(1003)
  obs.split <- group_initial_split(obs.ls[[i]], group=year)
  obs.train <- training(obs.split)
  obs.test <- testing(obs.split)
  
  
  
  # . prep ------------------------------------------------------------------
  
  responses <- c(alert="alert", tl="tl")[1]
  prep.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude))
  d.y <- list(train=map(prep.ls, ~bake(.x, obs.train)),
              test=map(prep.ls, ~bake(.x, obs.test)))
  prepPCA.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude, TRUE))
  dPCA.y <- list(train=map(prepPCA.ls, ~bake(.x, obs.train)),
                 test=map(prepPCA.ls, ~bake(.x, obs.test)))
  # saveRDS(prep.ls, glue("data/0_init/data_prep_{y}_{covSet}.rds"))
  # saveRDS(d.y, glue("data/0_init/data_baked_{y}_{covSet}.rds"))
  d.split <- dPCA.split <- map(responses, ~obs.split)
  for(r in responses) {
    d.split[[r]]$data <- d.split[[r]]$data %>% select(obsid) %>%
      left_join(bind_rows(d.y$train[[r]], d.y$test[[r]]))
    dPCA.split[[r]]$data <- dPCA.split[[r]]$data %>% select(obsid) %>%
      left_join(bind_rows(dPCA.y$train[[r]], dPCA.y$test[[r]]))
  }
  
  covs <- filter_corr_covs(all_covs, d.y, covs_exclude)
  covsPCA <- names(dPCA.y$train[[1]] %>% select(starts_with("PC")))
  if(any(table(d.y$train$tl$tl)==0)) {
    responses <- grep("tl", responses, invert=T, value=T)
  }
  
  # formulas
  smooths <- list(b=glue("b{c(covs$main, covs$interact)}"),
                  p=glue("p{c(covs$main, covs$interact)}"))
  smooths_PCA <- list(b=glue("b{covsPCA}"),
                      p=glue("p{covsPCA}"))
  form.ls <- map(
    responses,
    ~list(HBL=make_HB_formula(.x, c(covs$main, covs$interact)),
          HBL_PCA=make_HB_formula(.x, covsPCA),
          HBN=make_HB_formula(.x, c(covs$main, covs$interact), sTerms=smooths),
          HBN_PCA=make_HB_formula(.x, covsPCA, sTerms=smooths_PCA),
          ML=formula(glue("{.x} ~ {paste(unlist(covs), collapse='+')}")),
          ML_PCA=formula(glue("{.x} ~ {paste(covsPCA, collapse='+')}")),
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
                                 HBL_PCA=make_HB_priors(priStr, "HBL", .x, covsPCA, PCA=T),
                                 HBN=make_HB_priors(priStr, "HBN", .x, covs),
                                 HBN_PCA=make_HB_priors(priStr, "HBN", .x, covsPCA, PCA=T)))
  
  # tuning controls
  n_tuneVal <- list(Ridge=1e3,
                    ENet=1e3, 
                    MARS=1e3,
                    RF=1e2, 
                    NN=1e2,
                    Boost=1e2)
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
    set.seed(3001)
    folds <- vfold_cv(d.y$train[[r]], strata=r)
    set.seed(3001)
    foldsPCA <- vfold_cv(dPCA.y$train[[r]], strata=r)
    fit_model("Ridge", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("Ridge", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    fit_model("MARS", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("MARS", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    fit_model("RF", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("RF", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    fit_model("NN", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("NN", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    fit_model("Boost", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("Boost", r, form.ls, dPCA.y$train, foldsPCA, n_tuneVal, fit.dir, y, "_PCA")
    fit_model("HBL", r, form.ls, d.y$train, HB.i, priors, fit.dir, y)
    fit_model("HBL", r, form.ls, dPCA.y$train, HB.i, priors, fit.dir, y, "_PCA")
    fit_model("HBN", r, form.ls, d.y$train, HB.i, priors, fit.dir, y)
    fit_model("HBN", r, form.ls, dPCA.y$train, HB.i, priors, fit.dir, y, "_PCA")
  }
  
  fit.ls <- map(responses, ~summarise_predictions(d.y$train, dPCA.y$train, .x, fit.dir, y_i.i))
  saveRDS(fit.ls, glue("{out.dir}/{y}_fit_ls.rds"))
  
  oos.ls <- map(responses, ~summarise_predictions(d.y$test, dPCA.y$test, .x, fit.dir, y_i.i))
  saveRDS(oos.ls, glue("{out.dir}/{y}_oos_ls.rds"))
  
  
  
  # . cross validation by year ----------------------------------------------
  
  # Only necessary for Bayesian models; ML models use CV in fitting process
  for(r in responses) {
    set.seed(3001)
    folds <- vfold_cv(d.y$train[[r]], strata=r)
    set.seed(3001)
    foldsPCA <- vfold_cv(dPCA.y$train[[r]], strata=r)
    run_Bayes_CV("HBL", foldsPCA, cv.dir, y, y_i.i, r, form.ls, HB.i, priors, PCA=T)
    run_Bayes_CV("HBL", folds, cv.dir, y, y_i.i, r, form.ls, HB.i, priors)
    run_Bayes_CV("HBN", foldsPCA, cv.dir, y, y_i.i, r, form.ls, HB.i, priors, PCA=T)
    run_Bayes_CV("HBN", folds, cv.dir, y, y_i.i, r, form.ls, HB.i, priors)
  }
  
  
  
  # . ensemble --------------------------------------------------------------
  
  fit.ls <- merge_pred_dfs(dirf("out/compiled", glue("{y}_fit_ls.rds"), recursive=T))
  oos.ls <- merge_pred_dfs(dirf("out/compiled", glue("{y}_oos_ls.rds"), recursive=T))
  cv.ls <- list(alert=full_join(
    merge_pred_dfs(dirf("out/model_fits", glue("{y}_.*_HBL_CV"), recursive=T), CV="HB"),
    merge_pred_dfs(dirf("out/model_fits", glue("{y}_.*_CV.rds"), recursive=T), CV="ML"), 
    by=c("y", "obsid"))
    )
  wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
  
  # cl <- makeCluster(cores_per_model)
  # registerDoParallel(cl)
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, wt.ls, .x, y_i.i, "wtmean"))
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, cv.ls, .x, y_i.i, "GLM_fit", ens.dir, 1e3))
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, cv.ls, .x, y_i.i, "RF_fit", ens.dir, 1e2))
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, cv.ls, .x, y_i.i, "HB_fit", ens.dir, 0.8))
  # stopCluster(cl)
  
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, wt.ls, .x, y_i.i, "wtmean"))
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, cv.ls, .x, y_i.i, "RF_oos", ens.dir))
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, cv.ls, .x, y_i.i, "HB_oos", ens.dir))
  
  saveRDS(fit.ls, glue("out/compiled/{y}_fit.rds"))
  saveRDS(oos.ls, glue("out/compiled/{y}_oos.rds"))
  saveRDS(wt.ls, glue("out/compiled/{y}_wt.rds"))
  
  
  
  # . null ------------------------------------------------------------------
  
  fit.ls <- readRDS(glue("out/compiled/{y}_fit.rds")) 
  oos.ls <- readRDS(glue("out/compiled/{y}_oos.rds"))
  
  null.ls <- map(responses, ~calc_null(fit.ls, .x))
  fit.ls <- map(null.ls, ~.x$obs.df)
  oos.ls <- map2(oos.ls, null.ls, 
                 ~left_join(.x %>% mutate(yday=yday(date)), .y$yday.df) %>% select(-yday)) %>%
    map2(., fit.ls, ~bind_cols(.x, .y %>% select(contains("nullGrand")) %>% slice_head(n=1)))
  
  saveRDS(fit.ls, glue("out/compiled/{y}_fit.rds"))
  saveRDS(oos.ls, glue("out/compiled/{y}_oos.rds"))
  saveRDS(map(null.ls, ~.x$yday.df), glue("out/compiled/{y}_null.rds"))
  
}

closeAllConnections()

