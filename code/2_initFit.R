# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model fits



# setup -------------------------------------------------------------------
pkgs <- c("tidyverse", "lubridate", "glue", "recipes", "brms", "caret", 
          "randomForest", "RRF", "glmnet")
lapply(pkgs, library, character.only=T)
library(doParallel)
library(foreach)
source("code/00_fn.R")

cores_per_model <- 4
n_spp_parallel <- 7
test_startDate <- "2021-01-01"
fit.dir <- "out/model_fits/"
cv.dir <- "out/model_fits/cv/"

sp_i <- read_csv("data/sp_i.csv") %>% arrange(abbr)

col_metadata <- c("obsid", "sp", "date", "siteid", "lon", "lat")
col_resp <- c("lnN", "tl", "alert")
col_cmems <- readRDS("data/cmems_vars.rds")
col_wrf <- readRDS("data/wrf_vars.rds")

all_covs <- list(
  date=c("ydayCos", "ydaySin"),
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
  nonHB=c("ydaySinXydayCos", "lonz", "latz", "lonzXlatz")
)

obs.ls <- readRDS("data/0_init/data_full_allSpp.rds") %>%
  select(all_of(col_metadata), all_of(col_resp), 
         "alert1", "alert2", any_of(unname(unlist(all_covs)))) %>%
  mutate(across(starts_with("alert"), 
                ~factor(as.numeric(.x!="0_none"), labels=c("A0", "A1"))),
         across(starts_with("tl"), 
                ~factor(.x, ordered=T, 
                        levels=c("0_green", "1_yellow", "2_orange", "3_red"),
                        labels=c("TL0", "TL1", "TL2", "TL3")))) %>%
  mutate(year=year(date),
         yday=yday(date)) %>%
  group_by(sp) %>%
  group_split()
train.ls <- map(obs.ls, ~.x %>% filter(date < test_startDate))
test.ls <- map(obs.ls, ~.x %>% filter(date >= test_startDate))



# runs by taxon -----------------------------------------------------------

registerDoParallel(n_spp_parallel)
foreach(s=seq_along(train.ls),
        .export=c("all_covs", "train.ls", "test.ls", "fit.dir", "sp_i")
) %dopar% {
  
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  sp <- train.ls[[s]]$sp[1]
  sp_i.i <- filter(sp_i, abbr==sp)
  
  

# . prep ------------------------------------------------------------------

  # TODO: remove ifelse for final version
  reprep <- F
  responses <- c(alert="alert", tl="tl", lnN="lnN")
  if(reprep) {
    prep.ls <- map(responses, ~prep_recipe(train.ls[[s]], .x))
    d.sp <- list(train=map(prep.ls, ~bake(.x, train.ls[[s]])),
                    test=map(prep.ls, ~bake(.x, test.ls[[s]])))
    saveRDS(prep.ls, glue("data/0_init/data_prep_{sp}.rds"))
    saveRDS(d.sp, glue("data/0_init/data_baked_{sp}.rds"))
  } else {
    d.sp <- readRDS(glue("data/0_init/data_baked_{sp}.rds"))
  }
  
  covs <- filter_corr_covs(all_covs, d.sp, test_run=T)
  
  # formulas
  smooths <- list(b=glue("b{c(covs$main, covs$interact)}"),
                  p=glue("p{c(covs$main, covs$interact)}"))
  form.ls <- map(
    responses,
    ~list(HBL=make_HB_formula(.x, c(covs$main, covs$interact, covs$spacetime)),
          HBN=make_HB_formula(.x, c(covs$main, covs$interact), sTerms=smooths),
          ML=formula(glue("{.x} ~ {paste(unlist(covs), collapse='+')}"))
    )
  )
  
  # priors
  priStr <- switch(2,
                   "1"=list(r1=0.5, r2=2, hs1=0.5, hs2=0.6, b=0.75, de=0.3, i="1-loose"),
                   "2"=list(r1=0.3, r2=2, hs1=3, hs2=0.2, b=0.2, de=0.1, i="2-medium"),
                   "3"=list(r1=0.1, r2=2, hs1=5, hs2=0.5, b=0.5, de=0.05, i="3-tight"),
                   "4"=list(r1=0.3, r2=2, hs1=3, hs2=0.2, b=0.2, de=0.1, i="best")
  ) 
  priors <- map(responses, ~list(HBL=make_HB_priors(.x, priStr, "HBL"),
                                 HBN=make_HB_priors(.x, priStr, "HBN")))
  
  # tuning controls
  folds <- map(d.sp$train, createFoldsByYear)
  ctrl <- map(folds, ~trainControl("cv", classProbs=T, number=length(.x$i.in),
                                   index=.x$i.in, indexOut=.x$i.out))
  grids <- list(
    ENet=expand.grid(alpha=seq(0, 1, length.out=51),
                      lambda=2^(seq(-15,-1,length.out=50))),
    RRF=expand.grid(mtry=seq(1, min(3, length(unlist(covs))/10), by=1),
                   coefReg=seq(0.05, 0.8, by=0.05))
  )
  HB.i <- list(
    iter=1000,
    warmup=500,
    refresh=1,
    chains=cores_per_model,
    cores=cores_per_model,
    ctrl=list(adapt_delta=0.8, max_treedepth=10),
    prior_i=priStr$i
  )
  
  
  
# . train: fit models -----------------------------------------------------
  
  walk(responses, ~fit_model("ENet", .x, form.ls, d.sp$train, ctrl, grids, fit.dir, sp))
  walk(responses, ~fit_model("RRF", .x, form.ls, d.sp$train, ctrl, grids, fit.dir, sp))
  walk(responses, ~fit_model("HBL", .x, form.ls, d.sp$train, HB.i, priors, fit.dir, sp))
  walk(responses, ~fit_model("HBN", .x, form.ls, d.sp$train, HB.i, priors, fit.dir, sp))

  fit.ls <- map(responses, ~summarise_predictions(d.sp, "train", .x, fit.dir, sp_i.i))
  fit.ls$alert <- fit.ls$alert %>%
    full_join(fit.ls$tl %>% select(sp, obsid, ends_with("_A1"))) %>%
    full_join(fit.ls$lnN %>% select(sp, obsid, ends_with("_A1")))
  saveRDS(fit.ls, glue("out/{sp}_fit_ls.rds"))

  oos.ls <- map(responses, ~summarise_predictions(d.sp, "test", .x, fit.dir, sp_i.i))
  oos.ls$alert <- oos.ls$alert %>%
    full_join(oos.ls$tl %>% select(sp, obsid, ends_with("_A1"))) %>%
    full_join(oos.ls$lnN %>% select(sp, obsid, ends_with("_A1")))
  saveRDS(oos.ls, glue("out/{sp}_oos_ls.rds"))

  
  

# . cross validation by year ----------------------------------------------

  yrCV <- unique(d.sp$train$alert$year)
  for(k in 1:length(yrCV)) {
    yr <- yrCV[k]
    y_ <- paste0("_", yr)
    
    d.cv <- list(train=map(d.sp$train, ~.x %>% filter(year!=yr)),
                 test=map(d.sp$train, ~.x %>% filter(year==yr)))
    
    # tuning controls
    folds <- map(d.cv$train, createFoldsByYear)
    ctrl <- map(folds, ~trainControl("cv", classProbs=T, number=length(.x$i.in),
                                     index=.x$i.in, indexOut=.x$i.out))
    
    # fit models
    walk(responses, ~fit_model("ENet", .x, form.ls, d.cv$train, ctrl, grids, cv.dir, sp, y_))
    walk(responses, ~fit_model("RRF", .x, form.ls, d.cv$train, ctrl, grids, cv.dir, sp, y_))
    walk(responses, ~fit_model("HBL", .x, form.ls, d.cv$train, HB.i, priors, cv.dir, sp, y_))
    walk(responses, ~fit_model("HBN", .x, form.ls, d.cv$train, HB.i, priors, cv.dir, sp, y_))
    
    # predict
    cv.ls <- map(responses, ~summarise_predictions(d.cv, "test", .x, cv.dir, sp_i.i, y_))
    cv.ls$alert <- cv.ls$alert %>%
      full_join(cv.ls$tl %>% select(sp, obsid, ends_with("_A1"))) %>%
      full_join(cv.ls$lnN %>% select(sp, obsid, ends_with("_A1")))
    saveRDS(cv.ls, glue("{cv.dir}/{sp}_CV{y_}.rds"))
  }
  
  
  

# . ensemble --------------------------------------------------------------
  
  fit.ls <- readRDS(glue("out/{sp}_fit_ls.rds"))
  oos.ls <- readRDS(glue("out/{sp}_oos_ls.rds"))
  cv.ls <- map(dirf(cv.dir, glue("{sp}_CV")), readRDS)
  cv.ls <- list(alert=map_dfr(cv.ls, ~.x$alert),
                tl=map_dfr(cv.ls, ~.x$tl) %>% select(-ends_with("_A1")),
                lnN=map_dfr(cv.ls, ~.x$lnN) %>% select(-ends_with("_A1")))
  wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, wt.ls, .x, sp_i.i))
  fit.ls$alert <- fit.ls$alert %>%
    mutate(ens_tl_A1=fit.ls$tl$ens_tl_A1,
           ens_lnN_A1=fit.ls$lnN$ens_lnN_A1,
           ensAlt_lnN_A1=fit.ls$lnN$ensAlt_lnN_A1)
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, wt.ls, .x, sp_i.i))
  oos.ls$alert <- oos.ls$alert %>%
    mutate(ens_tl_A1=oos.ls$tl$ens_tl_A1,
           ens_lnN_A1=oos.ls$lnN$ens_lnN_A1,
           ensAlt_lnN_A1=oos.ls$lnN$ensAlt_lnN_A1)
  saveRDS(fit.ls, glue("out/{sp}_fit_ls.rds"))
  saveRDS(oos.ls, glue("out/{sp}_oos_ls.rds"))
  saveRDS(wt.ls, glue("out/{sp}_wt_ls.rds"))
  
  
  
  

# . null ------------------------------------------------------------------

  fit.ls <- readRDS(glue("out/{sp}_fit_ls.rds")) 
  oos.ls <- readRDS(glue("out/{sp}_oos_ls.rds"))
  
  null.ls <- map(responses, ~calc_null(fit.ls, .x))
  fit.ls <- map(null.ls, ~.x$obs.df)
  oos.ls <- map2(oos.ls, null.ls, 
                 ~left_join(.x %>% mutate(yday=yday(date)), .y$yday.df) %>% select(-yday))
  
  saveRDS(fit.ls, glue("out/{sp}_fit_ls.rds"))
  saveRDS(oos.ls, glue("out/{sp}_oos_ls.rds"))
  saveRDS(map(null.ls, ~.x$yday.df), glue("out/{sp}_null_ls.rds"))

  
  
  
  
  
}

closeAllConnections()




# Load dataset
# Identify predictors
# For each species:
# - scale variables that need it, storing the mean/sd 
# - make directional wind variables
# - make predictors factors that need it
# - split into testing/training
# - write_csv()
# Generate model formulas
# - lnN, alert
# Generate model priors
# Fit each model with training dataset
# CV by year
# Predict testing dataset
# Calculate ensemble
# Store predictions














# To Try ------------------------------------------------------------------

# Responses: lnN, alertBinary
# Interactions: cosDay:sinDay, lat:lon, fetch (?)
# Autocorrelation lags: blnN1 ~ lnDayLag1 + s(cosDay,sinDay) + s(lat,lon)
# Models: lnN
# - Null: 4wk-historic domain
# - Null: 4wk-historic site
# - ENetGLM (hurdle log-normal)
# - HB-int (hurdle log-normal)
# - HB-smooth (hurdle log-normal)
# - RRF
# Models: alertBinary
# - Null: 4wk-historic domain
# - Null: 4wk-historic site
# - ENetGLM
# - HB-int (ordinal)
# - HB-int (logistic)
# - HB-smooth (ordinal)
# - HB-smooth (logistic)
# - RRF
# - SVM
# - XGBoost
# - ANN

