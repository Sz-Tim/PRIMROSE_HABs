# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model fits



# setup -------------------------------------------------------------------
pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", "nnet", "randomForest", "glmnet", 
          "brms", "bayesian", "doParallel", "foreach")
lapply(pkgs, library, character.only=T)
source("code/00_fn.R")

covSet <- c("test", "full", "noDt", 
            "noDtDelta", "noInt", "noDtDeltaInt")[1]
cores_per_model <- 3
n_spp_parallel <- 6
fit.dir <- glue("out/model_fits/{covSet}/")
cv.dir <- glue("{fit.dir}/cv/")
out.dir <- glue("out/{covSet}/")
dir.create(cv.dir, recursive=T, showWarnings=F)
dir.create(out.dir, recursive=T, showWarnings=F)

y_i <- bind_rows(read_csv("data/i_hab.csv") %>% arrange(abbr) %>% mutate(type="hab"),
                 read_csv("data/i_tox.csv") %>% arrange(abbr) %>% mutate(type="tox"))

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

covs_exclude <- switch(covSet,
                       full="NA",
                       test="Xfetch|Dt|Delta|Avg",
                       noDt="Dt",
                       noDtDelta="Dt|Delta",
                       noInt="Xfetch",
                       noDtDeltaInt="Xfetch|Dt|Delta")

obs.ls <- map_dfr(dirf("data/0_init", "data_.*_all.rds"), readRDS) %>%
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
                  "fit.dir", "out.dir", "cv.dir", "y_i")
) %dopar% {
  
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  y <- obs.ls[[i]]$y[1]
  y_i.i <- filter(y_i, abbr==y)
  set.seed(123)
  obs.split <- group_initial_split(obs.ls[[i]], group=year)
  obs.train <- training(obs.split)
  obs.test <- testing(obs.split)
  
  

# . prep ------------------------------------------------------------------

  responses <- c(alert="alert", tl="tl")
  prep.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude))
  d.y <- list(train=map(prep.ls, ~bake(.x, obs.train)),
              test=map(prep.ls, ~bake(.x, obs.test)))
  # saveRDS(prep.ls, glue("data/0_init/data_prep_{y}_{covSet}.rds"))
  # saveRDS(d.y, glue("data/0_init/data_baked_{y}_{covSet}.rds"))
  d.split <- map(responses, ~obs.split)
  for(r in responses) {
    d.split[[r]]$data <- d.split[[r]]$data %>% select(obsid) %>%
      left_join(bind_rows(d.y$train[[r]], d.y$test[[r]]))
  }
  
  covs <- filter_corr_covs(all_covs, d.y, covs_exclude)
  if(any(table(d.y$train$tl$tl)==0)) {
    responses <- grep("tl", responses, invert=T, value=T)
  }
  
  # formulas
  smooths <- list(b=glue("b{c(covs$main, covs$interact)}"),
                  p=glue("p{c(covs$main, covs$interact)}"))
  form.ls <- map(
    responses,
    ~list(HBL=make_HB_formula(.x, c(covs$main, covs$interact)),
          HBN=make_HB_formula(.x, c(covs$main, covs$interact), sTerms=smooths),
          ML=formula(glue("{.x} ~ {paste(unlist(covs), collapse='+')}")),
          HB_vars=formula(glue("{.x} ~ yday + lon + lat + siteid +",
                               "{paste(unlist(covs), collapse='+')}"))
    )
  )
  
  # priors
  priStr <- switch(1,
                   "1"=list(r1=0.2, r2=2, hs1=0.5, hs2=0.6, b=0.75, de=0.3, i=1),
                   "2"=list(r1=0.1, r2=2, hs1=3, hs2=0.2, b=0.2, de=0.1, i=2),
                   "3"=list(r1=0.1, r2=1, hs1=5, hs2=0.5, b=0.5, de=0.05, i=3)
  ) 
  priors <- map(responses, ~list(HBL=make_HB_priors(priStr, "HBL", .x, covs),
                                 HBN=make_HB_priors(priStr, "HBN", .x, covs)))
  
  # tuning controls
  n_tuneVal <- 10
  HB.i <- list(
    iter=10,
    warmup=5,
    refresh=1,
    chains=cores_per_model,
    cores=cores_per_model,
    ctrl=list(adapt_delta=0.8, max_treedepth=10),
    prior_i=priStr$i
  )
  
  
  
# . train: fit models -----------------------------------------------------
  
  for(r in responses) {
    set.seed(123)
    folds <- vfold_cv(d.y$train[[r]], strata=r)
    fit_model("ENet", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("RF", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("NN", r, form.ls, d.y$train, folds, n_tuneVal, fit.dir, y)
    fit_model("HBL", r, form.ls, d.y$train, HB.i, priors, fit.dir, y)
    fit_model("HBN", r, form.ls, d.y$train, HB.i, priors, fit.dir, y)
  }
  
  fit.ls <- map(responses, ~summarise_predictions(d.y$train, .x, fit.dir, y_i.i))
  fit.ls$alert <- full_join(
    fit.ls$alert, 
    fit.ls[-which(responses=="alert")] %>%
      map(~.x %>% select(y, obsid, siteid, date, ends_with("_A1"))) %>%
      reduce(full_join)
  )
  saveRDS(fit.ls, glue("{out.dir}/{y}_fit_ls.rds"))
  
  oos.ls <- map(responses, ~summarise_predictions(d.y$test, .x, fit.dir, y_i.i))
  oos.ls$alert <- full_join(
    oos.ls$alert, 
    oos.ls[-which(responses=="alert")] %>%
      map(~.x %>% select(y, obsid, siteid, date, ends_with("_A1"))) %>%
      reduce(full_join)
  )
  saveRDS(oos.ls, glue("{out.dir}/{y}_oos_ls.rds"))
  
  

# . cross validation by year ----------------------------------------------

  yrCV <- unique(d.y$train$alert$year)
  for(k in 1:length(yrCV)) {
    yr <- yrCV[k]
    yr_ <- paste0("_", yr)
    
    d.cv <- list(train=map(d.y$train, ~.x %>% filter(year!=yr)),
                 test=map(d.y$train, ~.x %>% filter(year==yr)))
    
    # tuning controls
    folds <- map(d.cv$train, createFoldsByYear)
    ctrl <- map(folds, ~trainControl("cv", classProbs=T, number=length(.x$i.in),
                                     index=.x$i.in, indexOut=.x$i.out))
    
    for(r in responses) {
      set.seed(123)
      folds <- vfold_cv(d.cv$train[[r]])
      fit_model("ENet", r, form.ls, d.cv$train, folds, n_tuneVal, cv.dir, y, yr_)
      fit_model("RF", r, form.ls, d.cv$train, folds, n_tuneVal, cv.dir, y, yr_)
      fit_model("NN", r, form.ls, d.cv$train, folds, n_tuneVal, cv.dir, y, yr_)
      fit_model("HBL", r, form.ls, d.cv$train, HB.i, priors, cv.dir, y, yr_)
      fit_model("HBN", r, form.ls, d.cv$train, HB.i, priors, cv.dir, y, yr_)
    }
    
    # predict
    cv.ls <- map(responses, ~summarise_predictions(d.cv$test, .x, cv.dir, y_i.i, yr_))
    cv.ls$alert <- oos.ls %>%
      map(~.x %>% select(y, obsid, siteid, date, ends_with("_A1"))) %>%
      reduce(full_join)
    saveRDS(cv.ls, glue("{cv.dir}/{y}_CV{yr_}.rds"))
  }
  
  

# . ensemble --------------------------------------------------------------
  
  fit.ls <- readRDS(glue("{out.dir}/{y}_fit_ls.rds"))
  oos.ls <- readRDS(glue("{out.dir}/{y}_oos_ls.rds"))
  cv.ls <- map(dirf(cv.dir, glue("{y}_CV")), readRDS)
  cv.ls <- list(alert=map_dfr(cv.ls, ~.x$alert),
                tl=map_dfr(cv.ls, ~.x$tl) %>% select(-ends_with("_A1")),
                lnN=map_dfr(cv.ls, ~.x$lnN) %>% select(-ends_with("_A1")))
  wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
  
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, wt.ls, .x, y_i.i))
  fit.ls$alert <- full_join(
    fit.ls[[1]], 
    fit.ls[-1] %>%
      map(~.x %>% select(y, obsid, siteid, date, matches("ens_.*_A1"))) %>%
      reduce(full_join)
    )
  
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, wt.ls, .x, y_i.i))
  oos.ls$alert <- full_join(
    oos.ls[[1]], 
    oos.ls[-1] %>%
      map(~.x %>% select(y, obsid, siteid, date, matches("ens_.*_A1"))) %>%
      reduce(full_join)
  )
  
  saveRDS(fit.ls, glue("{out.dir}/{y}_fit_ls.rds"))
  saveRDS(oos.ls, glue("{out.dir}/{y}_oos_ls.rds"))
  saveRDS(wt.ls, glue("{out.dir}/{y}_wt_ls.rds"))
  
  

# . null ------------------------------------------------------------------

  fit.ls <- readRDS(glue("{out.dir}/{y}_fit_ls.rds")) 
  oos.ls <- readRDS(glue("{out.dir}/{y}_oos_ls.rds"))
  
  null.ls <- map(responses, ~calc_null(fit.ls, .x))
  fit.ls <- map(null.ls, ~.x$obs.df)
  oos.ls <- map2(oos.ls, null.ls, 
                 ~left_join(.x %>% mutate(yday=yday(date)), .y$yday.df) %>% select(-yday)) %>%
    map2(., fit.ls, ~bind_cols(.x, .y %>% select(contains("nullGrand")) %>% slice_head(n=1)))
  
  saveRDS(fit.ls, glue("{out.dir}/{y}_fit_ls.rds"))
  saveRDS(oos.ls, glue("{out.dir}/{y}_oos_ls.rds"))
  saveRDS(map(null.ls, ~.x$yday.df), glue("{out.dir}/{y}_null_ls.rds"))

  }

closeAllConnections()

