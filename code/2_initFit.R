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
    paste("VWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X")
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
                        labels=c("TL0green", "TL1yellow", "TL2orange", "TL3red")))) %>%
  mutate(year=year(date),
         ydayCos=cos(yday(date)),
         ydaySin=sin(yday(date))) %>%
  group_by(sp) %>%
  group_split()
train.ls <- map(obs.ls, ~.x %>% filter(date < test_startDate))
test.ls <- map(obs.ls, ~.x %>% filter(date >= test_startDate))




# runs by taxon -----------------------------------------------------------

registerDoParallel(n_spp_parallel)
foreach(s=seq_along(train.ls),
        .export=c("all_covs", "obs.ls", "train.ls", "test.ls")
) %dopar% {
  
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  sp <- train.ls[[s]]$sp[1]
  
  

# . prep ------------------------------------------------------------------

  # TODO: remove ifelse for final version
  reprep <- F
  responses <- c(lnN="lnN", tl="tl", alert="alert")
  if(reprep) {
    prep.ls <- map(responses, ~prep_recipe(train.ls[[s]], .x))
    data.sp <- list(train=map(prep.ls, ~bake(.x, train.ls[[s]])),
                    test=map(prep.ls, ~bake(.x, test.ls[[s]])))
    saveRDS(prep.ls, glue("data/0_init/data_prep_{sp}.rds"))
    saveRDS(data.sp, glue("data/0_init/data_baked_{sp}.rds"))
  } else {
    prep.ls <- readRDS(glue("data/0_init/data_prep_{sp}.rds"))
    data.sp <- readRDS(glue("data/0_init/data_baked_{sp}.rds"))
  }
  
  covs <- filter_corr_covs(all_covs, data.sp)
  
  # formulas
  smooths <- list(b=glue("b{c(covs$main, covs$interact)}"),
                  p=glue("p{c(covs$main, covs$interact)}"))
  form.ls <- map(
    responses,
    ~list(HBL=make_HB_formula(.x, c(all_covs$date, covs$main, covs$interact)),
          HBN=make_HB_formula(.x, c(covs$main, covs$interact), sTerms=smooths),
          GLM=formula(glue("{.x} ~ {paste(unlist(covs), collapse='+')}")),
          ML=formula(glue("{.x} ~ {paste(unlist(covs), collapse='+')}"))
    )
  )
  
  # priors
  priStr <- switch(2,
                   "1"=list(hs1=0.5, hs2=0.6, b=0.75, de=0.3, i="1-loose"),
                   "2"=list(hs1=3, hs2=0.2, b=0.2, de=0.1, i="2-medium"),
                   "3"=list(hs1=5, hs2=0.5, b=0.5, de=0.05, i="3-tight"),
                   "4"=list(hs1=3, hs2=0.2, b=0.2, de=0.1, i="best")
  ) 
  priors <- list(
    HBL=c(prior_string(glue("horseshoe({priStr$hs1}, par_ratio={priStr$hs2})"), class="b"),
          prior(normal(0, 1), class="Intercept"),
          prior(normal(0, 0.1), class="sd")),
    HBN=c(
      prior_string("normal(0,1)", class="b", nlpar="bIntercept"),
      prior_string("normal(0,.5)", class="sd", nlpar="bIntercept", lb=0),
      prior_string(glue("double_exponential(0,{priStr$de})"), class="sds", nlpar="bIntercept", lb=0),
      map(c(covs$main, covs$interact), 
          ~c(prior_string(glue("beta({priStr$b},1)"), nlpar=paste0("p", .x), lb=0, ub=1),
             prior_string("normal(0,1)", class="b", nlpar=paste0("b", .x)),
             prior_string("normal(0,.5)", class="sd", nlpar=paste0("b", .x), lb=0),
             prior_string("normal(0,.5)", class="sd", nlpar=paste0("p", .x), lb=0),
             prior_string(glue("double_exponential(0,{priStr$de})"), class="sds", nlpar=paste0("b", .x), lb=0))) %>%
        do.call('c', .)
    )
  )
  
  # Tuning controls
  folds <- map(data.sp$train, createFoldsByYear)
  ctrl <- map(folds, ~trainControl("cv", index=.x$i.in, indexOut=.x$i.out, classProbs=T))
  grids <- list(
    elast=expand.grid(alpha=seq(0, 1, length.out=26),
                      lambda=2^(seq(-15,-1,length.out=25))),
    rf=expand.grid(mtry=seq(1, min(6, length(unlist(covs))/10), by=1),
                   coefReg=seq(0.1, 0.8, by=0.1))
  )
  nTreeRF <- 20
  HB.i <- list(
    iter=100,
    warmup=50,
    refresh=1,
    chains=cores_per_model,
    cores=cores_per_model,
    ctrl=list(adapt_delta=0.95, max_treedepth=20)
  )
  
  
# . train: fit models -----------------------------------------------------
  
  # Elastic Net GLM
  elast.lnN <- train(form.ls$lnN$GLM, data=data.sp$train$lnN, 
                     method="glmnet", trControl=ctrl$lnN, tuneGrid=grids$elast)
  saveRDS(elast.lnN, glue("{fit.dir}/{sp}_lnN_elast.rds"))
  elast.tl <- train(form.ls$tl$GLM, data=data.sp$train$tl, 
                    method="glmnet", trControl=ctrl$tl, tuneGrid=grids$elast)
  saveRDS(elast.tl, glue("{fit.dir}/{sp}_tl_elast.rds"))
  elast.alert <- train(as.formula(form.ls$alert$GLM), data=data.sp$train$alert, 
                       method="glmnet", trControl=ctrl$alert, tuneGrid=grids$elast)
  saveRDS(elast.alert, glue("{fit.dir}/{sp}_alert_elast.rds"))
  
  # Random Forest
  rf.lnN <- train(form.ls$lnN$ML, data=data.sp$train$lnN, 
                    method="RRFglobal", trControl=ctrl$lnN, tuneGrid=rf.grid, ntree=nTreeRF)
  saveRDS(rf.lnN, glue("{fit.dir}/{sp}_lnN_rf.rds"))
  rf.tl <- train(form.ls$tl$ML, data=data.sp$train$tl, 
                 method="RRFglobal", trControl=ctrl$tl, tuneGrid=rf.grid, ntree=nTreeRF)
  saveRDS(rf.tl, glue("{fit.dir}/{sp}_tl_rf.rds"))
  rf.alert <- train(form.ls$alert$ML, data=data.sp$train$alert, 
                    method="RRFglobal", trControl=ctrl$alert, tuneGrid=rf.grid, ntree=nTreeRF)
  saveRDS(rf.lnN, glue("{fit.dir}/{sp}_alert_rf.rds"))
  
  # HB: linear
  HBL.lnN <- brm(form.ls$lnN$HBL, data=data.sp$train$lnN, family=hurdle_lognormal(), 
                 prior=priors$HBL, iter=HB.i$iter, warmup=HB.i$warmup, refresh=HB.i$refresh, 
                 init=0, control=HB.i$ctrl, chains=HB.i$chains, cores=HB.i$cores,
                 file=glue("{fit.dir}/{sp}_lnN_HBL.rds"))
  HBL.tl <- brm(form.ls$tl$HBL, data=data.sp$train$tl, family=cumulative(), 
                 prior=priors$HBL, iter=HB.i$iter, warmup=HB.i$warmup, refresh=HB.i$refresh, 
                 init=0, control=HB.i$ctrl, chains=HB.i$chains, cores=HB.i$cores,
                 file=glue("{fit.dir}/{sp}_tl_HBL.rds"))
  HBL.alert <- brm(form.ls$alert$HBL, data=data.sp$train$alert, family=bernoulli(), 
                 prior=priors$HBL, iter=HB.i$iter, warmup=HB.i$warmup, refresh=HB.i$refresh, 
                 init=0, control=HB.i$ctrl, chains=HB.i$chains, cores=HB.i$cores,
                 file=glue("{fit.dir}/{sp}_alert_HBL.rds"))
  
  # HB: nonlinear
  HBN.lnN <- brm(form.ls$lnN$HBN, data=data.sp$train$lnN, family=hurdle_lognormal(), 
                 prior=priors$HBN, iter=HB.i$iter, warmup=HB.i$warmup, refresh=HB.i$refresh, 
                 init=0, control=HB.i$ctrl, chains=HB.i$chains, cores=HB.i$cores,
                 file=glue("{fit.dir}/{sp}_lnN_HBN.rds"))
  HBN.tl <- brm(form.ls$tl$HBN, data=data.sp$train$tl, family=cumulative(), 
                prior=priors$HBN, iter=HB.i$iter, warmup=HB.i$warmup, refresh=HB.i$refresh, 
                init=0, control=HB.i$ctrl, chains=HB.i$chains, cores=HB.i$cores,
                file=glue("{fit.dir}/{sp}_tl_HBN.rds"))
  HBN.alert <- brm(form.ls$alert$HBN, data=data.sp$train$alert, family=bernoulli(), 
                   prior=priors$HBN, iter=HB.i$iter, warmup=HB.i$warmup, refresh=HB.i$refresh, 
                   init=0, control=HB.i$ctrl, chains=HB.i$chains, cores=HB.i$cores,
                   file=glue("{fit.dir}/{sp}_alert_HBN.rds"))
  
  
  # TEMP CHECKING
  varImp(elast.lnN)
  varImp(rf.lnN)
  par(mfrow=c(1,2))
  walk(data.sp, ~plot(predict(elast.lnN, .x$lnN), .x$lnN$lnN, col=rgb(0,0,0,0.2)))
  walk(data.sp, ~plot(predict(rf.lnN, .x$lnN), .x$lnN$lnN, col=rgb(0,0,0,0.2)))
  map(data.sp, ~yardstick::rmse_vec(.x$lnN$lnN, predict(elast.lnN, .x$lnN)))
  map(data.sp, ~yardstick::rmse_vec(.x$lnN$lnN, predict(rf.lnN, .x$lnN)))
  
  varImp(elast.tl)
  varImp(rf.tl)
  par(mfrow=c(4,2))
  walk(data.sp, ~plot(.x$tl$tl, predict(elast.tl, .x$tl, type="prob")[,1], main="GLM", ylim=c(0,1)))
  walk(data.sp, ~plot(.x$tl$tl, predict(elast.tl, .x$tl, type="prob")[,2], main="GLM", ylim=c(0,1)))
  walk(data.sp, ~plot(.x$tl$tl, predict(elast.tl, .x$tl, type="prob")[,3], main="GLM", ylim=c(0,1)))
  walk(data.sp, ~plot(.x$tl$tl, predict(elast.tl, .x$tl, type="prob")[,4], main="GLM", ylim=c(0,1)))
  walk(data.sp, ~plot(.x$tl$tl, predict(rf.tl, .x$tl, type="prob")[,1], main="RF", ylim=c(0,1)))
  walk(data.sp, ~plot(.x$tl$tl, predict(rf.tl, .x$tl, type="prob")[,2], main="RF", ylim=c(0,1)))
  walk(data.sp, ~plot(.x$tl$tl, predict(rf.tl, .x$tl, type="prob")[,3], main="RF", ylim=c(0,1)))
  walk(data.sp, ~plot(.x$tl$tl, predict(rf.tl, .x$tl, type="prob")[,4], main="RF", ylim=c(0,1)))
  map(data.sp, ~tibble(truth=.x$tl$tl) %>%
        bind_cols(predict(elast.tl, .x$tl, type="prob")) %>%
        yardstick::roc_auc(TL0green, TL1yellow, TL2orange, TL3red, truth=truth))
  map(data.sp, ~tibble(truth=.x$tl$tl) %>%
        bind_cols(predict(rf.tl, .x$tl, type="prob")) %>%
        yardstick::roc_auc(TL0green, TL1yellow, TL2orange, TL3red, truth=truth))
  map(data.sp, ~tibble(truth=.x$tl$tl) %>%
        bind_cols(predict(elast.tl, .x$tl, type="prob")) %>%
        yardstick::roc_curve(TL0green, TL1yellow, TL2orange, TL3red, truth=truth) %>%
        ggplot(aes(1-specificity, sensitivity, colour=.level)) + geom_line() + ggtitle("GLM") +
        scale_colour_manual(values=c("green3", "yellow3", "orange", "red3")))
  map(data.sp, ~tibble(truth=.x$tl$tl) %>%
        bind_cols(predict(rf.tl, .x$tl, type="prob")) %>%
        yardstick::roc_curve(TL0green, TL1yellow, TL2orange, TL3red, truth=truth) %>%
        ggplot(aes(1-specificity, sensitivity, colour=.level)) + geom_line() + ggtitle("RF") +
        scale_colour_manual(values=c("green3", "yellow3", "orange", "red3")))
  
  par(mfrow=c(1,2))
  varImp(elast.alert)
  varImp(rf.alert)
  walk(data.sp, ~plot(.x$alert$alert, predict(elast.alert, .x$alert, type="prob")[,2]))
  walk(data.sp, ~plot(.x$alert$alert, predict(rf.alert, .x$alert, type="prob")[,2]))
  map(data.sp, ~yardstick::roc_auc_vec(.x$alert$alert, 
                                       predict(elast.alert, .x$alert, type="prob")[,2],
                                       event_level="second"))
  map(data.sp, ~yardstick::roc_auc_vec(.x$alert$alert, 
                                       predict(rf.alert, .x$alert, type="prob")[,2],
                                       event_level="second"))
  
  
  
  # Random Forest
  
  
  
  

# . train: summarise ------------------------------------------------------
  
  fit.lnN <- data.sp$train$lnN %>%
    select(sp, obsid, siteid, date, lnN) %>%
    mutate(
      glm_lnN=predict(elast.lnN),
      rf_lnN=predict(rf.lnN),
      HBL_lnN,
      HBN_lnN
    )
  fit.tl <- data.sp$train$tl %>%
    select(sp, obsid, siteid, date, tl) %>%
    mutate(
      glm_tl=predict(elast.tl),
      rf_tl=predict(rf.tl),
      HBL_tl,
      HBN_tl
    )
  fit.alert <- data.sp$train$alert %>%
    select(sp, obsid, siteid, date, alert) %>%
    mutate(
      glm_alert=predict(elast.alert),
      rf_alert=predict(rf.alert),
      HBL_alert,
      HBN_alert
    )
  

# . test ------------------------------------------------------------------


  

# . cross-validate --------------------------------------------------------


  

# . ensemble --------------------------------------------------------------

      
  
}






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












lnN.recipe <- recipe(lnN ~ fetch + phycWk + phycDt + lnNAvg1 + ydaySin + ydayCos + lat + lon, obs.ls[[1]]) %>%
  step_normalize(any_of(col_cmems))
lnN.wf <- workflow() %>%
  add_model(rand_forest("regression") %>% set_engine("ranger")) %>%
  add_recipe(lnN.recipe)
test <- lnN.wf %>%
  fit(obs.ls[[1]])
test.pred <- predict(test, obs.ls[[2]])
plot(obs.ls[[2]]$lnN, test.pred$.pred)





obs.df %>%
  filter(date < test_startDate) %>%
  arrange(sp, date, siteid) %>%
  saveRDS("data/0_init/data_allSpp_train.rds")
obs.df %>%
  filter(date >= test_startDate) %>%
  arrange(sp, date, siteid) %>%
  saveRDS("data/0_init/data_allSpp_test.rds")







# To Try ------------------------------------------------------------------

# Responses: lnN, alertBinary
# Interactions: cosDay:sinDay, lat:lon, fetch (?)
# Autocorrelation lags: blnN1 ~ lnDayLag1 + s(cosDay,sinDay) + s(lat,lon)
# Models: lnN
# - Null: 4wk-historic domain
# - Null: 4wk-historic site
# - ElastGLM (hurdle log-normal)
# - HB-int (hurdle log-normal)
# - HB-smooth (hurdle log-normal)
# - RF
# Models: alertBinary
# - Null: 4wk-historic domain
# - Null: 4wk-historic site
# - ElastGLM
# - HB-int (ordinal)
# - HB-int (logistic)
# - HB-smooth (ordinal)
# - HB-smooth (logistic)
# - RF
# - SVM
# - XGBoost
# - ANN

