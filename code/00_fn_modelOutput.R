# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Functions: model output





#' Summarise and aggregate predictions for a set of models
#'
#' @param d.y 
#' @param set 
#' @param resp 
#' @param fit.dir 
#' @param y_i.i 
#' @param suffix 
#'
#' @return
#' @export
#'
#' @examples
summarise_predictions <- function(d.y, dPCA.y, resp, fit.dir, y_i.i, suffix=NULL) {
  library(tidyverse); library(glue); library(tidymodels)
  fits.f <- dirf(fit.dir, glue("{y_i.i$abbr[1]}_{resp}.*{ifelse(is.null(suffix),'',suffix)}"))
  names(fits.f) <- str_split_fixed(str_split_fixed(fits.f, glue("{resp}_"), 2)[,2], "_|\\.", 2)[,1]
  fits <- map(fits.f, readRDS)
  
  fits.dPCA <- grep("PCA", fits.f)
  fits.d <- grep("PCA", fits.f, invert=T)
  
  if(length(fits.dPCA) > 0) {
    preds.dPCA <- imap_dfc(fits[fits.dPCA], ~get_predictions(.x, .y, resp, dPCA.y, y_i.i)) |>
      rename_with(~glue("PCA.{.x}"))
  } else {
    preds.dPCA <- NULL
  }
  
  if(length(fits.d) > 0) {
    preds.d <- imap_dfc(fits[fits.d], ~get_predictions(.x, .y, resp, d.y, y_i.i))
  } else {
    preds.d <- NULL
  }
  
  out.df <- d.y[[resp]] |>
    select(y, obsid, siteid, date, prevAlert, {{resp}}) 
  if(!is.null(preds.d)) {
    out.df <- out.df |> bind_cols(preds.d)
  }
  if(!is.null(preds.dPCA)) {
    out.df <- out.df |> bind_cols(preds.dPCA)
  }
  return(out.df)
}





#' Generate predictions for a model
#'
#' @param fit 
#' @param mod 
#' @param resp 
#' @param d.df 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
get_predictions <- function(fit, mod, resp, d.df, y_i.i) {
  library(tidyverse); library(glue); library(tidymodels); library(brms)
  
  if(grepl("HB", mod)) {
    pred <- parsnip::extract_fit_engine(fit) |>
      posterior_epred(d.df[[resp]], allow_new_levels=T) |>
      summarise_post_preds(resp, y_i.i)
  } else {
    pred_type <- ifelse(resp=="lnN", "raw", "prob")
    preds <- predict(fit, d.df[[resp]], pred_type)
    pred <- summarise_ML_preds(preds, resp, y_i.i)
  }
  pred.df <- as_tibble(pred) |> 
    rename_with(~glue("{mod}_{resp}_{.x}"))
  return(pred.df)
}




#' Calculate mean predictions for ML models
#'
#' @param preds 
#' @param resp 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
summarise_ML_preds <- function(preds, resp, y_i.i) {
  if(resp=="alert") {
    pred <- cbind(A1=preds$.pred_A1)
  }
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1)) + 1
    pred <- cbind(preds, 
                  rowSums(preds[,thresh:4]))
    colnames(pred) <- c(paste0("TL", 0:3), "A1")
  }
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    pred <- cbind(lnN=preds,
                  A1=preds>=thresh)
  }
  return(pred)
}




#' Calculate mean predictions for HB models
#'
#' @param post 
#' @param resp 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
summarise_post_preds <- function(post, resp, y_i.i) {
  if(resp=="alert") {
    pred <- cbind(A1=colMeans(post))
  }
  
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1)) + 1
    pred <- cbind(apply(post, 2:3, mean),
                  colMeans(apply(post[,,(thresh):4], 1:2, sum)))
    colnames(pred) <- c(paste0("TL", 0:3), "A1")
  }
  
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    pred <- cbind(lnN=colMeans(post),
                  A1=colMeans(post>=thresh))
  }
  return(pred)
} 









#' Merge summarised predictions across all models into single dataframes
#'
#' @param files vector of full file path
#'
#' @return
#' @export
#'
#' @examples
merge_pred_dfs <- function(files, CV=NULL) {
  f.df <- tibble(f=files, 
                 covSet=str_split(files, "/") |> 
                   map_chr(~grep("^[0-9]", .x, value=T) |> 
                             str_split_fixed("-", 2) |> 
                             magrittr::extract(,1)))
  if(is.null(CV)) {
    map(1:nrow(f.df), 
        ~readRDS(f.df$f[.x]) |> 
          lapply(function(x) {x |> mutate(covSet=paste0("d", f.df$covSet[.x], "."))})) |>
      list_transpose() |>
      map_depth(2, ~.x |> 
                  pivot_longer(ends_with("_A1"), names_to="model", values_to="prA1") |>
                  mutate(model=paste0(covSet, model)) |>
                  select(-covSet) |>
                  pivot_wider(names_from="model", values_from="prA1")) |>
      map(~reduce(.x, full_join)) 
  } else if(CV=="HB") {
    map_dfr(1:nrow(f.df), 
            ~readRDS(f.df$f[.x]) |> mutate(covSet=paste0("d", f.df$covSet[.x], "."))) |>
      pivot_longer(ends_with("A1"), names_to="model", values_to="prA1") |>
      select(-prevAlert) |>
      na.omit() |>
      mutate(model=paste0(covSet, model)) |>
      select(-covSet) |>
      pivot_wider(names_from="model", values_from="prA1")
  } else if(CV=="ML") {
    map(1:nrow(f.df), 
        ~readRDS(f.df$f[.x]) |> 
          mutate(covSet=paste0("d", f.df$covSet[.x], ".")) |>
          pivot_longer(ends_with("A1"), names_to="model", values_to="prA1") |>
          mutate(model=paste0(covSet, model)) |>
          select(-covSet) |>
          pivot_wider(names_from="model", values_from="prA1")) |>
      reduce(full_join)
  } else {
    stop("CV must be 'HB', 'ML', or NULL")
  }
}








#' Calculate model weights based on mean log loss
#'
#' @param cv.df 
#' @param resp 
#' @param wt.penalty 
#'
#' @return
#' @export
#'
#' @examples
calc_LL_wts <- function(cv.df, resp, wt.penalty=2) {
  library(yardstick)
  if(resp=="alert") {
    wt.df <- cv.df |>
      pivot_longer(ends_with("_A1"), names_to="model", values_to="pr") |>
      group_by(model) |>
      average_precision(pr, truth=alert, event_level="second") |>
      mutate(.estimate=log(.estimate))
  }
  if(resp=="tl") {
    wt.df <- cv.df |>
      pivot_longer(contains("_tl_TL"), names_to="ID", values_to="pr") |>
      mutate(model=str_split_fixed(ID, "_TL", 2)[,1],
             predCat=str_split_fixed(ID, "_tl_", 2)[,2]) |>
      select(-ID) |>
      pivot_wider(names_from=predCat, values_from=pr) |>
      group_by(model) |>
      mn_log_loss(TL0, TL1, TL2, TL3, truth=tl)
  }
  if(resp=="lnN") {
    wt.df <- cv.df |>
      pivot_longer(ends_with("_lnN_lnN"), names_to="model", values_to="pr") |>
      group_by(model) |>
      rmse(truth=lnN, estimate=pr)
  }
  return(wt.df |>
           ungroup() |>
           mutate(wt=(1/.estimate^wt.penalty)/sum(1/.estimate^wt.penalty)))
}





#' Calculate ensemble model predictions
#'
#' @param out.ls 
#' @param wt.ls 
#' @param resp 
#' @param y_i.i 
#'
#' @return
#' @export
#'
#' @examples
calc_ensemble <- function(out.ls, wt.ls, resp, y_i.i, method="wtmean", out.path=NULL, opt=NULL) {
  library(tidyverse); library(tidymodels)
  if(grepl("RF|GLM|HB", method) & resp != "alert") {
    stop("Models only implemented for alert")
  }
  if(resp=="alert") {
    if(method=="wtmean") {
      out <- left_join(
        out.ls[[resp]] |> 
          pivot_longer(ends_with("_A1"), names_to="model", values_to="pr"),
        wt.ls[[resp]]
      ) |>
        mutate(pr_logit=brms::logit_scaled(pr, lb=-0.01, ub=1.01)) |>
        group_by(obsid) |>
        summarise(ens_alert_A1=sum(pr*wt, na.rm=T),
                  ensLogitMn_alert_A1=brms::inv_logit_scaled(
                    sum(pr_logit*wt, na.rm=T))) |>
        ungroup()
    } else if(grepl("[GLM|RF|HB]_fit", method)) {
      avg_prec2 <- metric_tweak("avg_prec2", average_precision, event_level="second")
      folds <- vfold_cv(wt.ls[[resp]], strata="alert")
      ens_rec <- recipe(alert~., data=wt.ls[[resp]]) |>
        update_role(y, obsid, siteid, date, new_role="ID") |>
        step_logit(ends_with("_A1"), offset=0.01) |>
        step_normalize(all_predictors())
      ens_rec2 <- recipe(alert~., data=wt.ls[[resp]]) |>
        update_role(y, obsid, siteid, date, new_role="ID") 
      
      if(method=="GLM_fit") {
        size <- ifelse(is.null(opt), 1e3, opt)
        GLM_spec <- logistic_reg(penalty=tune(), mixture=0) |>
          set_engine("glmnet", lower.limits=0) |> set_mode("classification")
        
        GLM_wf <- workflow() |>
          add_model(GLM_spec) |>
          add_recipe(ens_rec)
        GLM_tune <- GLM_wf |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(GLM_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        GLM_out <- GLM_wf |>
          finalize_workflow(select_best(GLM_tune, "avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(GLM_out, glue("{out.path}/{y_i.i$abbr}_EnsGLM.rds"))
        
        GLM_wf2 <- workflow() |>
          add_model(GLM_spec) |>
          add_recipe(ens_rec2)
        GLM_tune2 <- GLM_wf2 |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(GLM_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        GLM_out2 <- GLM_wf2 |>
          finalize_workflow(select_best(GLM_tune2, "avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(GLM_out2, glue("{out.path}/{y_i.i$abbr}_EnsGLM2.rds"))
      }
      
      if(method=="RF_fit") {
        size <- ifelse(is.null(opt), 1e2, opt)
        RF_spec <- rand_forest(trees=tune(), 
                               min_n=tune()) |>
          set_engine("randomForest") |> set_mode("classification")
        
        RF_wf <- workflow() |>
          add_model(RF_spec) |>
          add_recipe(ens_rec)
        RF_tune <- RF_wf |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(RF_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        RF_out <- RF_wf |>
          finalize_workflow(select_best(RF_tune, "avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(RF_out, glue("{out.path}/{y_i.i$abbr}_EnsRF.rds"))
        
        RF_wf2 <- workflow() |>
          add_model(RF_spec) |>
          add_recipe(ens_rec2)
        RF_tune2 <- RF_wf2 |>
          tune_grid(resamples=folds,
                    grid=grid_latin_hypercube(extract_parameter_set_dials(RF_spec),
                                              size=size),
                    metrics=metric_set(avg_prec2))
        RF_out2 <- RF_wf2 |>
          finalize_workflow(select_best(RF_tune2, "avg_prec2")) |>
          fit(wt.ls[[resp]]) |>
          butcher()
        saveRDS(RF_out2, glue("{out.path}/{y_i.i$abbr}_EnsRF2.rds"))
      }
      
      if(method=="HB_fit") {
        mn <- ifelse(is.null(opt), 0.5, opt)
        library(brms)
        wf <- workflow() |>
          add_model(bayesian(mode="classification", engine="brms", 
                             family=bernoulli, 
                             prior=prior_string(glue("exponential({mn})"), class="b", lb=0),
                             init=0, 
                             iter=1500,
                             warmup=1000,
                             chains=3,
                             cores=3,
                             save_model=glue("{out.path}/{y_i.i$abbr}_EnsHB1.stan")),
                    formula=alert~.) |>
          add_recipe(ens_rec2)
        HB_out <- wf |>
          fit(data=wt.ls[[resp]]) 
        saveRDS(HB_out |> axe_env_bayesian() |> axe_env_bayesian(), 
                glue("{out.path}/{y_i.i$abbr}_EnsHB.rds"))
      }
    }
    if(grepl("GLM", method)) {
      GLM_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsGLM.rds"))
      GLM_out2 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsGLM2.rds"))
      out <- out.ls[[resp]] %>%
        mutate(ensGLM_alert_A1=predict(GLM_out, new_data=., type="prob")[[2]],
               ensGLM2_alert_A1=predict(GLM_out2, new_data=., type="prob")[[2]]) |>
        select(obsid, starts_with("ensGLM")) 
    }
    if(grepl("RF", method)) {
      RF_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsRF.rds"))
      RF_out2 <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsRF2.rds"))
      out <- out.ls[[resp]] %>%
        mutate(ensRF_alert_A1=predict(RF_out, new_data=., type="prob")[[2]],
               ensRF2_alert_A1=predict(RF_out2, new_data=., type="prob")[[2]]) |>
        select(obsid, starts_with("ensRF")) 
    }
    if(grepl("HB", method)) {
      HB_out <- readRDS(glue("{out.path}/{y_i.i$abbr}_EnsHB.rds"))
      pred <- parsnip::extract_fit_engine(HB_out) |>
        posterior_epred(out.ls[[resp]], allow_new_levels=T) |>
        summarise_post_preds(resp, y_i.i)
      out <- out.ls[[resp]] %>%
        mutate(ensHB_alert_A1=pred[,1]) |>
        select(obsid, starts_with("ensHB")) 
    }
  }
  if(resp=="tl") {
    thresh <- as.numeric(str_sub(y_i.i$tl_thresh, -1, -1))
    alert_cols <- paste0("ens_tl_TL", thresh:3)
    out <- left_join(
      out.ls[[resp]] |>
        select(-ends_with("_A1")) |>
        pivot_longer(contains("_tl_TL"), names_to="ID", values_to="pr") |>
        mutate(model=str_split_fixed(ID, "_TL", 2)[,1],
               predCat=str_split_fixed(ID, "_tl_", 2)[,2]) |>
        select(-ID),
      wt.ls[[resp]]) |>
      group_by(obsid, predCat) |>
      summarise(ens_tl=sum(pr*wt, na.rm=T)) |>
      ungroup() |>
      pivot_wider(names_from=predCat, values_from=ens_tl, names_prefix="ens_tl_") |>
      mutate(ens_tl_A1=rowSums(pick(all_of(alert_cols))))
  }
  if(resp=="lnN") {
    thresh <- log1p(y_i.i$N_thresh)
    out <- left_join(
      out.ls[[resp]] |>
        select(-ends_with("_A1")) |>
        pivot_longer(ends_with("_lnN_lnN"), names_to="model", values_to="pr"),
      wt.ls[[resp]]) |>
      group_by(obsid) |>
      summarise(ens_lnN=sum(pr*wt, na.rm=T),
                ens_lnN_A1=sum((pr >= thresh)*wt, na.rm=T)) |>
      ungroup()
  }
  return(left_join(out.ls[[resp]], out))
}






#' Identify rows in y within thresh days of x 
#'
#' @param x 
#' @param y 
#' @param thresh 
#'
#' @return
#' @export
#'
#' @examples
get_rows_by_yday <- function(x, y, thresh) {
  xy_diff <- x - y
  comp.mx <- cbind(abs(xy_diff),
                   abs(xy_diff - 365),
                   abs(xy_diff + 365))
  ids <- which(apply(comp.mx, 1, function(j) any(j <= thresh)))
  return(ids)
}






#' Calculate null model predictions
#'
#' @param obs.ls 
#' @param resp 
#'
#' @return
#' @export
#'
#' @examples
calc_null <- function(obs.ls, resp) {
  library(tidyverse)
  obs.df <- obs.ls[[resp]] |>
    mutate(yday=yday(date)) |>
    select(yday, {{resp}})
  if(resp=="alert") {
    temp.df <- tibble(yday=1:365,
                      null4wk_alert_A1=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_alert_A1[i] <- mean(obs.df$alert[obs.ids] == "A1")
    }
    out.df <- obs.ls[[resp]] |> 
      mutate(yday=yday(date),
             nullGrand_alert_A1=mean(alert=="A1"))
  }
  if(resp=="tl") {
    temp.df <- tibble(yday=1:365,
                      null4wk_tl_TL0=NA_real_,
                      null4wk_tl_TL1=NA_real_,
                      null4wk_tl_TL2=NA_real_,
                      null4wk_tl_TL3=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_tl_TL0[i] <- mean(obs.df$tl[obs.ids] == "TL0")
      temp.df$null4wk_tl_TL1[i] <- mean(obs.df$tl[obs.ids] == "TL1")
      temp.df$null4wk_tl_TL2[i] <- mean(obs.df$tl[obs.ids] == "TL2")
      temp.df$null4wk_tl_TL3[i] <- mean(obs.df$tl[obs.ids] == "TL3")
    }
    out.df <- obs.ls[[resp]] |> 
      mutate(yday=yday(date),
             nullGrand_tl_TL0=mean(tl=="TL0"),
             nullGrand_tl_TL1=mean(tl=="TL1"),
             nullGrand_tl_TL2=mean(tl=="TL2"),
             nullGrand_tl_TL3=mean(tl=="TL3"))
  }
  if(resp=="lnN") {
    temp.df <- tibble(yday=1:365,
                      null4wk_lnN_lnN=NA_real_)
    for(i in 1:nrow(temp.df)) {
      obs.ids <- get_rows_by_yday(i, obs.df$yday, 14)
      temp.df$null4wk_lnN_lnN[i] <- mean(obs.df$lnN[obs.ids])
    }
    out.df <- obs.ls[[resp]] |> 
      mutate(yday=yday(date),
             nullGrand_lnN_lnN=mean(lnN))
  }
  return(list(yday.df=temp.df,
              obs.df=left_join(out.df, temp.df) |> select(-yday)))
}









#' Calculate J, F-1, precision, and recall for a dataframe
#' 
#' This is a function to pass to future_map to avoid capturing large objects
#' in the environment which makes it extraordinarily slow
#'
#' @param dat 
#'
#' @return
#' @export
#'
#' @examples
get_metrics <- function(dat) {
  tibble(F1=f_meas(dat, alert, pred, event_level="second")$.estimate,
         F2=f_meas(dat, alert, pred, beta=2, event_level="second")$.estimate,
         precision=precision(dat, alert, pred, event_level="second")$.estimate,
         recall=recall(dat, alert, pred, event_level="second")$.estimate,
         accuracy=accuracy(dat, alert, pred)$.estimate,)
}









#' Compute F-1, J-index, precision, and recall across probability thresholds
#'
#' @param L.df 
#' @param prSteps 
#'
#' @return
#' @export
#'
#' @examples
compute_thresholds <- function(L.df, prMin=0, prMax=1, prSteps=0.1, byPrevAlert=F, cores=1) {
  library(tidyverse); library(yardstick); library(furrr)
  pred.df <- map_dfr(seq(prMin, prMax, by=prSteps), 
                     ~L.df |> mutate(thresh=.x)) |>
    mutate(pred=if_else(prA1 < thresh, "A0", "A1") |> factor(levels=c("A0", "A1")))
  col_to_drop <- c("obsid", "siteid", "date", "prA1") 
  if(!byPrevAlert) {
    col_to_drop <- c(col_to_drop, "prevAlert")
  }
  if(cores > 1) {
    plan(multisession, workers=cores)
  } else {
    plan(sequential)
  }
  metric.df <- pred.df |>
    select(-all_of(col_to_drop)) |>
    nest(dat=c(alert, pred)) |>
    ungroup() |>
    mutate(metrics=future_map(dat, get_metrics)) |>
    select(-dat) |>
    unnest(metrics)
  return(metric.df)
}









#' Find minimum possible PR-AUC
#' 
#' Following Boyd et al 2012
#'
#' @param df Dataframe
#' @param ... Unquoted grouping variables
#'
#' @return
#' @export
#'
#' @examples
find_AUCPR_min <- function(df, ...) {
  df_aucpr_min <- df |> 
    filter(model==levels(model)[1]) |>
    group_by(...) |>
    summarise(prop=mean(alert=="A1")) |>
    mutate(AUCPR_min=1+((1-prop)*log(1-prop))/prop) |>
    ungroup() |>
    select(-prop)
  return(left_join(df, df_aucpr_min))
}




