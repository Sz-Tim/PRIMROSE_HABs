# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model fits



# setup -------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(glue)
library(recipes)
library(brms)
library(caret)
library(randomForest)
library(xgboost)
library(e1071)
library(glmnet)
library(doParallel)
source("code/00_fn.R")

cores_per_model <- 4
n_spp_parallel <- 7
test_startDate <- "2021-01-01"
fit.dir <- "out/model_fits/"

sp_i <- read_csv("data/sp_i.csv") %>% arrange(abbr)

col_metadata <- c("obsid", "sp", "date", "siteid", "lon", "lat")
col_resp <- c("lnN", "tl", "alert")
col_cmems <- readRDS("data/cmems_includeVars.rds")
col_wrf <- c("U", "UDir", "Shortwave", "Precip", "sst")

all_covs <- list(
  date=c("ydayCos", "ydaySin"),
  main=c(
    "fetch", 
    "lnNWt1", "alert1A1", "lnNAvg1", "prAlertAvg1",
    "lnNWt2", "alert2A1", "lnNAvg2", "prAlertAvg2",
    col_cmems
  ),
  interact=c(),
  nonHB=c("ydaySinXydayCos", "lonz", "latz", "lonzXlatz")
)

obs.ls <- readRDS("data/0_init/data_allSpp_full.rds") %>%
  select(all_of(col_metadata), all_of(col_resp), any_of(unname(unlist(all_covs)))) %>%
  mutate(across(starts_with("alert"), ~factor(as.numeric(.x!="0_none"), labels=c("A0", "A1"))),
         across(starts_with("tl"), ~factor(.x, ordered=T))) %>%
  mutate(year=year(date),
         ydayCos=cos(yday(date)),
         ydaySin=sin(yday(date))) %>%
  group_by(sp) %>%
  group_split()
train.ls <- map(obs.ls, ~.x %>% filter(date < test_startDate))
test.ls <- map(obs.ls, ~.x %>% filter(date >= test_startDate))




# runs by taxon -----------------------------------------------------------
# parallelize with n_spp_parallel
for(s in seq_along(obs.ls)) {
  
  sp <- obs.ls[[s]]$sp[1]
  

# . prep ------------------------------------------------------------------

  prep.ls <- map(c(lnN="lnN", tl="tl", alert="alert"),
                 ~prep_recipe(train.ls[[s]], .x))
  data.sp <- list(train=map(prep.ls, ~bake(.x, train.ls[[s]])),
                  test=map(prep.ls, ~bake(.x, test.ls[[s]])))
  saveRDS(prep.ls, glue("data/0_init/data_prep_{sp}.rds"))
  saveRDS(data.sp, glue("data/0_init/data_baked_{sp}.rds"))
  
  covs <- filter_corr_covs(all_covs, data.sp)
  
  smooths <- list(b=glue("b{c(covs$main, covs$interact)}"),
                  p=glue("p{c(covs$main, covs$interact)}"))
  form.ls <- list(
    lnN=list(
      HB_lin=make_HB_formula("lnN", c(covs$date, covs$main, covs$interact)),
      HB_nl=make_HB_formula("lnN", c(covs$main, covs$interact), sTerms=smooths),
      GLM=formula(glue("lnN ~ {paste(unlist(covs), collapse='+')}")),
      ML=formula(glue("lnN ~ {paste(unlist(covs), collapse='+')}"))
    ),
    tl=list(
      HB_lin=make_HB_formula("tl", c(covs$date, covs$main, covs$interact)),
      HB_nl=make_HB_formula("tl", c(covs$main, covs$interact), sTerms=smooths)
    ),
    alert=list(
      HB_lin=make_HB_formula(resp="alert", c(covs$date, covs$main, covs$interact)),
      HB_nl=make_HB_formula("alert", c(covs$main, covs$interact), sTerms=smooths),
      GLM=formula(glue("alert ~ {paste(unlist(covs), collapse='+')}")),
      ML=formula(glue("alert ~ {paste(unlist(covs), collapse='+')}"))
    )
  )
  
  
# . train: fit models -----------------------------------------------------

  # Elastic Net GLM
  elast.ctrl <- trainControl("repeatedcv", number=6, repeats=3, classProbs=T)
  elast.grid <- expand.grid(alpha=seq(0, 1, length.out=26),
                            lambda=2^(seq(-15,-1,length.out=25)))
  elast.lnN <- train(form.ls$lnN$GLM, 
                     data=data.sp$train$lnN, 
                     method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
  saveRDS(elast.lnN, glue("{fit.dir}/{sp}_lnN_elast.rds"))
  elast.alert <- train(as.formula(form.ls$alert$GLM), 
                       data=data.sp$train$alert, 
                       method="glmnet", trControl=elast.ctrl, tuneGrid=elast.grid)
  saveRDS(elast.alert, glue("{fit.dir}/{sp}_alert_elast.rds"))
  
  # TEMP CHECKING
  plot(predict(elast.lnN, data.sp$train$lnN), data.sp$train$lnN$lnN, col=rgb(0,0,0,0.2))
  plot(predict(elast.alert, data.sp$train$alert, type="prob")[,1], data.sp$train$alert$alert, col=rgb(0,0,0,0.05))
  plot(predict(elast.lnN, data.sp$test$lnN), data.sp$test$lnN$lnN, col=rgb(0,0,0,0.2))
  plot(predict(elast.alert, data.sp$test$alert, type="prob")[,2], data.sp$test$alert$alert, col=rgb(0,0,0,0.05))
  yardstick::pr_auc_vec(data.sp$train$alert$alert, predict(elast.alert, data.sp$train$alert, type="prob")[,2])
  yardstick::pr_auc_vec(data.sp$test$alert$alert, predict(elast.alert, data.sp$test$alert, type="prob")[,2])
  # Random Forest
  
  
  
  

# . train: summarise ------------------------------------------------------

  fit.df <- full_join(
    train.df %>%
      mutate(glmElast_mnpr=predict(elast_, train.df_, type="prob")[,2]),
    tibble(obsid=c(filter(train.df, Nbloom1==0)$obsid, filter(train.df, Nbloom1==1)$obsid),
           glmElast_split_mnpr=c(
             predict(elast01, train.df01, type="prob")[,2],
             predict(elast11, train.df11, type="prob")[,2]))
  ) %>%
    mutate(covarSet=i.name,
           species=target)
  

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

