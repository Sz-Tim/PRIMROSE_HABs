# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Model periodic re-fits



# setup -------------------------------------------------------------------
pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", 
          "nnet", "randomForest", "glmnet", "xgboost", "earth",
          "brms", "bayesian", "doParallel", "foreach", "butcher")
lapply(pkgs, library, character.only=T)
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
n_spp_parallel <- 20
base.dir <- "out/" 




# Model predictions -------------------------------------------------------

# for(i in 1:nrow(covSet.df)) {
registerDoParallel(n_spp_parallel)
foreach(i=1:nrow(covSet.df)) %dopar% {
  
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  
  covSet <- covSet.df$f[i]
  d <- covSet.df$id[i]
  
  fit.dir <- glue("{base.dir}/model_fits/{covSet}/")
  cv.dir <- glue("{fit.dir}/cv/")
  ens.dir <- glue("{base.dir}/ensembles/")
  out.dir <- glue("{base.dir}/compiled/{covSet}/")
  dir.create(ens.dir, recursive=T, showWarnings=F)
  dir.create(cv.dir, recursive=T, showWarnings=F)
  dir.create(out.dir, recursive=T, showWarnings=F)
  
  y_i <- bind_rows(read_csv("data/i_hab.csv", show_col_types=F) |> 
                     arrange(abbr) |> mutate(type="hab"),
                   read_csv("data/i_tox.csv", show_col_types=F) |> 
                     arrange(abbr) |> mutate(type="tox")) |>
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
  prepPCA.ls <- map(responses, ~prep_recipe(obs.train, .x, covs_exclude, TRUE))
  dPCA.y <- list(train=map(prepPCA.ls, ~bake(.x, obs.train)),
                 test=map(prepPCA.ls, ~bake(.x, obs.test)))
  
  fit.ls <- map(responses, ~summarise_predictions(d.y$train, dPCA.y$train, .x, fit.dir, y_i.i))
  saveRDS(fit.ls, glue("{out.dir}/{y}_fit_ls.rds"))
  
  oos.ls <- map(responses, ~summarise_predictions(d.y$test, dPCA.y$test, .x, fit.dir, y_i.i))
  saveRDS(oos.ls, glue("{out.dir}/{y}_oos_ls.rds"))
  
}


closeAllConnections()




# Compile -----------------------------------------------------------------

n_spp_parallel <- 20
ens.dir <- glue("{base.dir}/ensembles/")
registerDoParallel(n_spp_parallel)
for(i in seq_along(y_resp)) {
# foreach(i=seq_along(y_resp),
#         .export=c("ens.dir", "y_i")
# ) %dopar% {
  
  lapply(pkgs, library, character.only=T)
  source("code/00_fn.R")
  y <- y_resp[i]
  y_i.i <- filter(y_i, abbr==y)
  set.seed(1003)
  
  responses <- c(alert="alert", tl="tl")[1]
  
  # . ensemble --------------------------------------------------------------
  
  fit.ls <- merge_pred_dfs(dirf("out/compiled", glue("{y}_fit_ls.rds"), recursive=T))
  oos.ls <- merge_pred_dfs(dirf("out/compiled", glue("{y}_oos_ls.rds"), recursive=T))
  cv.ls <- list(alert=full_join(
    merge_pred_dfs(dirf("out/model_fits", glue("{y}_.*_HB[L|N]_CV"), recursive=T), CV="HB"),
    merge_pred_dfs(dirf("out/model_fits", glue("{y}_.*_CV.rds"), recursive=T), CV="ML"),
    by=c("y", "obsid"))
  )
  cv.ls$alert <- cv.ls$alert |> na.omit()
  wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
  saveRDS(cv.ls, glue("out/compiled/{y}_cv.rds"))
  saveRDS(wt.ls, glue("out/compiled/{y}_wt.rds"))
  
  
  # cl <- makeCluster(cores_per_model)
  # registerDoParallel(cl)
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, wt.ls, .x, y_i.i, "wtmean"))
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, cv.ls, .x, y_i.i, "GLM_fit", ens.dir, 1e3))
  # fit.ls <- map(responses, ~calc_ensemble(fit.ls, cv.ls, .x, y_i.i, "RF_fit", ens.dir, 1e2))
  # stopCluster(cl)
  # fit.ls <- map(responses, ~calc_ensemble(fit.ls, cv.ls, .x, y_i.i, "HB_fit", ens.dir, 2))
  #
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, wt.ls, .x, y_i.i, "wtmean"))
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
  # oos.ls <- map(responses, ~calc_ensemble(oos.ls, cv.ls, .x, y_i.i, "RF_oos", ens.dir))
  # oos.ls <- map(responses, ~calc_ensemble(oos.ls, cv.ls, .x, y_i.i, "HB_oos", ens.dir))
  
  saveRDS(fit.ls, glue("out/compiled/{y}_fit.rds"))
  saveRDS(oos.ls, glue("out/compiled/{y}_oos.rds"))
  
  
  
  # . null ------------------------------------------------------------------
  
  fit.ls <- readRDS(glue("out/compiled/{y}_fit.rds"))
  oos.ls <- readRDS(glue("out/compiled/{y}_oos.rds"))
  
  null.ls <- map(responses, ~calc_null(fit.ls, .x))
  fit.ls <- map(null.ls, ~.x$obs.df)
  oos.ls <- map2(oos.ls, null.ls,
                 ~left_join(.x |> mutate(yday=yday(date)), .y$yday.df) |> select(-yday)) |>
    map2(.x=_, fit.ls, ~bind_cols(.x, .y |> select(contains("nullGrand")) |> slice_head(n=1)))
  
  saveRDS(fit.ls, glue("out/compiled/{y}_fit.rds"))
  saveRDS(oos.ls, glue("out/compiled/{y}_oos.rds"))
  saveRDS(map(null.ls, ~.x$yday.df), glue("out/compiled/{y}_null.rds"))
  
}

closeAllConnections()



