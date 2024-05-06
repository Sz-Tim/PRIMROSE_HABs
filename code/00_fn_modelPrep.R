# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Functions: model preparation




#' Get regex for excluded covariates given a covariate set
#'
#' @param covSet 
#'
#' @return
#' @export
#'
#' @examples
get_excluded_cov_regex <- function(covSet) {
  covSet.df <- expand_grid(Avg=c(0,1), 
                           Xf=c(0,1),
                           XN=c(0,1),
                           Del=c(0,1)) |>
    mutate(id=row_number(),
           f=glue("{id}-Avg{Avg}_Xf{Xf}_XN{XN}_Del{Del}"),
           exclude=glue("Dt|", 
                        "{ifelse(Avg==0, 'Avg|', 'NA|')}",
                        "{ifelse(Xf==0, 'Xfetch|', 'NA|')}",
                        "{ifelse(XN==0, 'lnNWt1X|', 'NA|')}",
                        "{ifelse(Del==0, 'Delta', 'NA')}"))
  filter(covSet.df, f==covSet)$exclude
}









#' Add Dt dummy variables
#' 
#' Should only be necessary for initial model fitting since there were some
#' where Dt variables were included but then excluded within the recipe. The
#' bake() function still expects them initially even though they are excluded.
#'
#' @param df 
#'
#' @return
#' @export
#'
#' @examples
add_dummy_columns <- function(df) {
  dummy_cols <- c("chlAvgDtDirE", "chlAvgDtDirN", "chlAvgDtDirS", "chlAvgDtDirW", "chlDt", 
                  "kdAvgDtDirE", "kdAvgDtDirN", "kdAvgDtDirS", "kdAvgDtDirW", "kdDt", 
                  "no3AvgDtDirE", "no3AvgDtDirN", "no3AvgDtDirS", "no3AvgDtDirW", "no3Dt", 
                  "o2AvgDtDirE", "o2AvgDtDirN", "o2AvgDtDirS", "o2AvgDtDirW", "o2Dt",
                  "phAvgDtDirE", "phAvgDtDirN", "phAvgDtDirS", "phAvgDtDirW", "phDt", 
                  "phycAvgDtDirE", "phycAvgDtDirN", "phycAvgDtDirS", "phycAvgDtDirW", "phycDt", 
                  "po4AvgDtDirE", "po4AvgDtDirN", "po4AvgDtDirS", "po4AvgDtDirW", "po4Dt", 
                  "ppAvgDtDirE", "ppAvgDtDirN", "ppAvgDtDirS", "ppAvgDtDirW", "ppDt", 
                  "PrecipAvgDtDirE", "PrecipAvgDtDirN", "PrecipAvgDtDirS", "PrecipAvgDtDirW", "PrecipDt",
                  "ShortwaveAvgDtDirE", "ShortwaveAvgDtDirN", "ShortwaveAvgDtDirS", "ShortwaveAvgDtDirW", "ShortwaveDt", 
                  "sstAvgDtDirE", "sstAvgDtDirN", "sstAvgDtDirS", "sstAvgDtDirW", "sstDt",
                  "UAvgDtDirE", "UAvgDtDirN", "UAvgDtDirS", "UAvgDtDirW", "UDt", 
                  "UVAvgDtDirE", "UVAvgDtDirN", "UVAvgDtDirS", "UVAvgDtDirW", "UVDt", 
                  "VAvgDtDirE", "VAvgDtDirN", "VAvgDtDirS", "VAvgDtDirW", "VDt")
  df |>
    bind_cols(map_dfr(dummy_cols, ~tibble(x=.x, y=0)) |> pivot_wider(names_from=x, values_from=y))
}








#' Create recipe and prepare using training data
#'
#' @param train.df 
#' @param response 
#' @param dimReduce 
#'
#' @return
#' @export
#'
#' @examples
prep_recipe <- function(train.df, response, covsExclude="NA", dimReduce=FALSE) {
  respExclude <- grep(response, c("lnN", "tl", "alert"), value=T, invert=T)
  pred_vars <- names(train.df)
  include_UVX <- !grepl("Xfetch", covsExclude)
  include_lnNX <- !grepl("lnNWt1X", covsExclude)
  rec <- recipe(train.df) |>
    step_mutate(prevAlert=alert1, role="ID") |>
    update_role(all_of(pred_vars), new_role="predictor") |>
    update_role(all_of(response), new_role="outcome") |>
    update_role(obsid, y, date, siteid, year, new_role="ID") |>
    update_role(lon, lat, new_role="RE") |>
    step_select(-any_of(respExclude)) |> 
    step_dummy(all_factor_predictors()) |>
    step_logit(starts_with("prAlert"), offset=0.01) |>
    step_logit(ends_with("A1"), offset=0.01) |>
    step_harmonic(yday, frequency=1, cycle_size=365, keep_original_cols=T) |>
    update_role(yday, new_role="RE") |>
    step_rename(ydayCos=yday_cos_1, ydaySin=yday_sin_1) |>
    step_mutate_at(lon, lat, fn=list(z=~.)) |>
    step_interact(term=~ydaySin:ydayCos, sep="X")
  if(include_UVX) {
    rec <- rec |>
      step_interact(terms=~UWk:fetch:matches("Dir[EW]"), sep="X") |>
      step_interact(terms=~VWk:fetch:matches("Dir[NS]"), sep="X") |>
      step_interact(terms=~UWk:fetch:matches("Dir[NS]"), sep="X") |>
      step_interact(terms=~VWk:fetch:matches("Dir[EW]"), sep="X")
  }
  if(include_lnNX) {
    rec <- rec |>
      step_interact(terms=~lnNWt1:all_predictors(), sep="X") 
  }
  rec <- rec |>
    step_interact(terms=~lon_z:lat_z, sep="X") |>
    step_YeoJohnson(all_predictors()) |>
    step_normalize(all_predictors()) |>
    step_rename_at(contains("_"), fn=~str_remove_all(.x, "_")) |>
    step_select(-matches(covsExclude))
  if(dimReduce) {
    rec <- rec |>
      step_pca(all_predictors(), threshold=0.95)
  }
  rec |>
    prep(training=train.df)
}






#' Filter list of covariates following recipe thinning
#'
#' @param all_covs 
#' @param data.y 
#' @param covsExclude 
#'
#' @return
#' @export
#'
#' @examples
filter_corr_covs <- function(all_covs, data.y, covsExclude="NA") {
  uncorr_covs <- unique(unlist(map(data.y, ~unlist(map(.x, names)))))
  if(any(grepl("^PC", uncorr_covs))) {
    covs_incl <- list(main=grep("^PC", uncorr_covs, value=T),
                      interact=NULL,
                      nonHB=NULL)
  } else {
    uncorr_covs <- unique(c(uncorr_covs, "ydayCos", "ydaySin", "lon", "lat"))
    covs_incl <- map(all_covs, ~.x[.x %in% uncorr_covs])
  }
  
  covs_incl <- map(covs_incl, ~grep(covsExclude, .x, value=T, invert=T))
  return(covs_incl)
}





#' Create formulas for Hierarchical Bayesian models
#'
#' @param resp 
#' @param covs 
#' @param sTerms 
#' @param splinesInt 
#' @param splinesCovs 
#'
#' @return
#' @export
#'
#' @examples
make_HB_formula <- function(resp, covs, sTerms=NULL, 
                            splinesInt="both", splinesCovs="time") {
  library(tidyverse); library(brms); library(glue)
  
  splines_int <- switch(splinesInt,
                        "time"="s(yday, bs=('cc'))",
                        "space"="s(lon, lat)",
                        "both"="t2(yday, lon, lat, bs=c('cc','ts'), d=c(1,2))")
  splines_cov <- switch(splinesCovs,
                        "time"="s(yday, bs=('cc'))",
                        "space"="s(lon, lat)",
                        "both"="t2(yday, lon, lat, bs=c('cc','ts'), d=c(1,2))")
  
  if(is.null(sTerms)) {
    if(resp=="lnN") {
      form <- bf(glue("{resp} ~ 1 + {paste(covs, collapse='+')}", 
                      "+ {splines_int}",
                      "+ (1 + {paste(covs, collapse='+')} | siteid)"),
                 glue("hu ~ 1 + {paste(covs, collapse='+')}", 
                      "+ {splines_int}",
                      "+ (1 + {paste(covs, collapse='+')} | siteid)"))  
    } else {
      form <- bf(glue("{resp} ~ 1 + {paste(covs, collapse='+')}", 
                      "+ {splines_int}",
                      "+ (1 + {paste(covs, collapse='+')} | siteid)"))
    }
  } else {
    flist <- c(
      map(sTerms$b, ~as.formula(glue("{.x} ~ {splines_cov} + (1|siteid)"))),
      map(sTerms$p, ~as.formula(glue("{.x} ~ 1 + (1|siteid)"))),
      map(1, ~glue("bIntercept ~ 1 + {splines_int} + (1|siteid)"))
    )
    if(resp=="lnN") {
      sTerms$pHu <- paste0(sTerms$p, "Hu")
      sTerms$bHu <- paste0(sTerms$b, "Hu")
      flist <- c(
        nlf(glue("hu ~ 1*bInterceptHu",
                 "+ {paste(sTerms$pHu, sTerms$bHu, covs, sep='*', collapse='+')}")),
        flist,
        map(sTerms$bHu, ~as.formula(glue("{.x} ~ {splines_cov} + (1|siteid)"))),
        map(sTerms$pHu, ~as.formula(glue("{.x} ~ 1 + (1|siteid)"))),
        map(1, ~glue("bInterceptHu ~ 1 + {splines_int} + (1|siteid)"))
      )
      form <- bf(glue("{resp} ~ 1*bIntercept",
                      "+ {paste(sTerms$p, sTerms$b, covs, sep='*', collapse='+')}"),
                 flist=flist, nl=T) 
    } else {
      form <- bf(glue("{resp} ~ 1*bIntercept",
                      "+ {paste(sTerms$p, sTerms$b, covs, sep='*', collapse='+')}"),
                 flist=flist, nl=T)
    }
  }
  return(form)
}







#' Create priors for each Hierarchical Bayesian model
#'
#' @param prior_i 
#' @param mod 
#' @param resp 
#' @param covs 
#'
#' @return
#' @export
#'
#' @examples
make_HB_priors <- function(prior_i, mod, resp, covs, PCA=F) {
  library(tidyverse); library(brms)
  if(mod=="HBL") {
    p <- c(prior(normal(0, 1), class="Intercept"),
           prior(normal(0, 0.1), class="sd"))
    if(resp=="tl") {
      p <- c(p,
             prior_string(glue("horseshoe({prior_i$hs1}, par_ratio={prior_i$hs2})"), class="b"))
    } else {
      p <- c(p,
             prior_string(glue("R2D2({prior_i$r1}, {prior_i$r2})"), class="b"))
    }
    if(resp=="lnN") {
      p <- c(p,
             prior_string(glue("R2D2({prior_i$r1}, {prior_i$r2})"), class="b", dpar="hu"))
    }
  }
  if(mod=="HBN") {
    if(PCA) {
      terms <- list(int="bIntercept",
                    cov=covs)
    } else {
      terms <- list(int="bIntercept",
                    cov=c(covs$main, covs$interact))
    }
    
    if(resp=="lnN") {
      terms <- map(terms, ~c(.x, paste0(.x, "Hu")))
    }
    p <- c(
      map(terms$int, 
          ~c(prior_string("normal(0,1)", class="b", nlpar=.x),
             prior_string("normal(0,.5)", class="sd", nlpar=.x, lb=0),
             prior_string("student_t(3,0,2.5)", class="sds", nlpar=.x, lb=0))) |>
        do.call('c', args=_),
      map(terms$cov, 
          ~c(prior_string(glue("beta({prior_i$b},1)"), 
                          nlpar=paste0("p", .x), lb=0, ub=1),
             prior_string("normal(0,1)", class="b", 
                          nlpar=paste0("b", .x)),
             prior_string("normal(0,.5)", class="sd", 
                          nlpar=paste0("b", .x), lb=0),
             prior_string("normal(0,.5)", class="sd", 
                          nlpar=paste0("p", .x), lb=0),
             prior_string(glue("double_exponential(0,{prior_i$de})"), class="sds", 
                          nlpar=paste0("b", .x), lb=0))) |>
        do.call('c', args=_)
    )
  }
  return(p)
}





#' Identify indexes for ML training
#'
#' @param data.df 
#'
#' @return
#' @export
#'
#' @examples
createFoldsByYear <- function(data.df) {
  folds_out <- data.df |> mutate(rowNum=row_number()) |> 
    group_by(year) |> group_split() |> map(~.x$rowNum)
  folds_out <- folds_out[-1]
  folds_in <- map(folds_out, ~(1:nrow(data.df))[-.x])
  return(list(i.in=folds_in, i.out=folds_out))
}







#' Wrapper to fit a model
#'
#' @param mod 
#' @param resp 
#' @param form.ls 
#' @param d.ls 
#' @param opts 
#' @param tunes 
#' @param out.dir 
#' @param y 
#' @param suffix 
#'
#' @return
#' @export
#'
#' @examples
fit_model <- function(mod, resp, form.ls, d.ls, opts, tunes, out.dir, y, suffix=NULL) {
  library(glue); library(tidymodels)
  dir.create(glue("{out.dir}/meta/"), showWarnings=F, recursive=T)
  dir.create(glue("{out.dir}/vi/"), showWarnings=F, recursive=T)
  PCA_run <- all(!is.null(suffix), grepl("PCA", suffix))
  # Fit ML models
  if(mod %in% c("Ridge", "ENet", "RF", "NN", "MARS", "Boost", "lgbm")) {
    fit_ID <- glue("{y}_{resp}_{mod}{ifelse(is.null(suffix),'',suffix)}")
    if(file.exists(glue("{out.dir}/{fit_ID}.rds"))) {
      cat("File already exists:", glue("{out.dir}/{fit_ID}.rds"), "\n")
      return()
    }
    mod.prefix <- ifelse(PCA_run, "PCA.", "")
    ML_form <- ifelse(PCA_run, "ML_PCA", "ML")
    ML_spec <- switch(mod,
                      Ridge=logistic_reg(penalty=tune(), 
                                         mixture=0) |>
                        set_engine("glmnet") |> set_mode("classification"),
                      ENet=logistic_reg(penalty=tune(), 
                                        mixture=tune()) |>
                        set_engine("glmnet") |> set_mode("classification"),
                      RF=rand_forest(trees=tune(), 
                                     min_n=tune()) |>
                        set_engine("randomForest") |> set_mode("classification"),
                      NN=mlp(hidden_units=tune(), 
                             penalty=tune(), 
                             epochs=tune()) |>
                        set_engine("nnet", maxNWts=1e4) |> set_mode("classification"),
                      MARS=mars(num_terms=tune(), 
                                prod_degree=tune()) |>
                        set_engine("earth") |> set_mode("classification"),
                      Boost=boost_tree(trees=tune(),
                                       tree_depth=tune(),
                                       min_n=tune(),
                                       learn_rate=tune(),
                                       loss_reduction=tune()) |>
                        set_engine("xgboost") |> set_mode("classification"),
                      lgbm=boost_tree(tree_depth=tune(),
                                      trees=tune(),
                                      learn_rate=tune(),
                                      mtry=tune(),
                                      min_n=tune(),
                                      loss_reduction=tune()) |>
                        set_engine("lightgbm") |> set_mode("classification")
    )
    avg_prec2 <- metric_tweak("avg_prec2", average_precision, event_level="second")
    pr_auc2 <- metric_tweak("pr_auc2", pr_auc, event_level="second")
    roc_auc2 <- metric_tweak("roc_auc2", roc_auc, event_level="second")
    wf <- workflow() |>
      add_model(ML_spec) |>
      add_formula(form.ls[[resp]][[ML_form]])
    set.seed(1003)
    if(mod=="lgbm") {
      out_tune <- wf |>
        tune_grid(resamples=opts, 
                  grid=grid_latin_hypercube(
                    extract_parameter_set_dials(ML_spec) |> 
                      update(mtry=finalize(mtry(), opts)), 
                    size=tunes[[mod]]),
                  metrics=metric_set(roc_auc2, pr_auc2, avg_prec2), 
                  control=control_grid(save_pred=T))  
    } else {
      out_tune <- wf |>
        tune_grid(resamples=opts, 
                  grid=grid_latin_hypercube(extract_parameter_set_dials(ML_spec), 
                                            size=tunes[[mod]]),
                  metrics=metric_set(roc_auc2, pr_auc2, avg_prec2), 
                  control=control_grid(save_pred=T))
      # saveRDS(out_tune |> butcher, glue("{out.dir}/meta/{fit_ID}_tune.rds")) 
    }
    best <- select_best(out_tune, metric="avg_prec2")
    out_tune |> 
      collect_predictions() |>
      filter(.config==best$.config) |>
      arrange(.row) |>
      mutate(obsid=d.ls[[resp]]$obsid[.row],
             y=y) |>
      select(y, obsid, .pred_A1) |>
      rename_with(~glue("{mod.prefix}{mod}_{resp}_A1"), .cols=".pred_A1") |>
      saveRDS(glue("{out.dir}/cv/{fit_ID}_CV.rds"))
    out <- wf |>
      finalize_workflow(best) |>
      fit(data=d.ls[[resp]]) 
    out |> 
      extract_fit_engine() |>
      vip::vi(scale=T) |>
      saveRDS(glue("{out.dir}/vi/{fit_ID}_vi.rds"))
    out <- out |>
      butcher()
  }
  
  # Fit Hierarchical Bayesian models
  if(mod %in% c("HBL", "HBN")) {
    library(brms)
    fit_ID <- glue("{y}_{resp}_{mod}{opts$prior_i}{ifelse(is.null(suffix),'',suffix)}")
    if(file.exists(glue("{out.dir}/{fit_ID}.rds"))) {
      cat("File already exists:", glue("{out.dir}/{fit_ID}.rds"), "\n")
      return()
    }
    HB.family <- switch(resp, 
                        lnN=hurdle_lognormal,
                        tl=cumulative,
                        alert=bernoulli)
    HB_form <- ifelse(PCA_run, paste0(mod, "_PCA"), paste0(mod))
    HB_form_dummy <- ifelse(PCA_run, "HB_vars_PCA", "HB_vars")
    wf <- workflow() |>
      add_model(bayesian(mode="classification", engine="brms", 
                         formula.override=bayesian_formula(form.ls[[resp]][[HB_form]]),
                         family=HB.family, 
                         prior=tunes[[resp]][[HB_form]],
                         init=0, 
                         iter=opts$iter,
                         warmup=opts$warmup,
                         control=opts$ctrl,
                         chains=opts$chains,
                         cores=opts$cores,
                         save_model=glue("{out.dir}/meta/{fit_ID}.stan")),
                formula=form.ls[[resp]]$HB_vars) |>
      add_recipe(recipe(d.ls[[resp]], formula=form.ls[[resp]][[HB_form_dummy]]))
    out <- wf |>
      fit(data=d.ls[[resp]])
    out |> 
      extract_fit_engine() |>
      fixef() |>
      as_tibble(rownames="Variable") |>
      saveRDS(glue("{out.dir}/vi/{fit_ID}_vi.rds"))
    out <- out |>
      axe_env_bayesian() |>
      axe_env_bayesian()
  }
  saveRDS(out, glue("{out.dir}/{fit_ID}.rds"))
  cat("Saved ", y, "_", resp, "_", mod, " as ", out.dir, "*", suffix, "\n", sep="")
}









#' Run cross-validation for Bayesian models
#'
#' @param mod 
#' @param folds 
#' @param cv.dir 
#' @param y 
#' @param y_i.i 
#' @param r 
#' @param form.ls 
#' @param HB.i 
#' @param priors 
#'
#' @return
#' @export
#'
#' @examples
run_Bayes_CV <- function(mod, folds, cv.dir, y, y_i.i, r, form.ls, HB.i, priors, PCA=F, reverse=F) {
  for(f in 1:nrow(folds)) {
    if(reverse) {
      f <- (nrow(folds):1)[f] 
    }
    f_ <- paste0("_f", str_pad(f, 2, side="left", pad="0"))
    f_ <- ifelse(PCA, paste0("_PCA", f_), f_)
    d.cv <- list(train=list(alert=training(folds$splits[[f]])),
                 test=list(alert=testing(folds$splits[[f]])))
    if(!file.exists(glue("{cv.dir}/{y}_{r}_{mod}_CV{f_}.rds"))) {
      fit_model(mod, r, form.ls, d.cv$train, HB.i, priors, cv.dir, y, f_)
    }
    if(file.exists(glue("{cv.dir}/{y}_{r}_{mod}1{f_}.rds"))) {
      summarise_predictions(d.cv$test, d.cv$test, r, cv.dir, y_i.i, glue("{mod}1{f_}")) |>
        saveRDS(glue("{cv.dir}/{y}_{r}_{mod}_CV{f_}.rds"))
      file.remove(glue("{cv.dir}/{y}_{r}_{mod}1{f_}.rds"))
    }
  }
}
