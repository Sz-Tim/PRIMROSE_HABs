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
base.dir <- "out/1_current/" 




# Model predictions -------------------------------------------------------

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
  
  obs.ls <- map_dfr(dirf("data/1_current", "data_.*_all.rds"), readRDS) |>
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
  
  fit.ls <- merge_pred_dfs(dirf("out/1_current/compiled", glue("{y}_fit_ls.rds"), recursive=T))
  oos.ls <- merge_pred_dfs(dirf("out/1_current/compiled", glue("{y}_oos_ls.rds"), recursive=T))
  cv.ls <- list(alert=full_join(
    merge_pred_dfs(dirf("out/1_current/model_fits", glue("{y}_.*_HB[L|N]_CV"), recursive=T), CV="HB"),
    merge_pred_dfs(dirf("out/1_current/model_fits", glue("{y}_.*_CV.rds"), recursive=T), CV="ML"),
    by=c("y", "obsid"))
  )
  cv.ls$alert <- cv.ls$alert |> na.omit()
  wt.ls <- imap(cv.ls, ~calc_LL_wts(.x, .y))
  saveRDS(cv.ls, glue("out/1_current/compiled/{y}_cv.rds"))
  saveRDS(wt.ls, glue("out/1_current/compiled/{y}_wt.rds"))
  
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, wt.ls, .x, y_i.i, "wtmean"))
  fit.ls <- map(responses, ~calc_ensemble(fit.ls, cv.ls, .x, y_i.i, "GLM_fit", ens.dir, 1e3))
  #
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, wt.ls, .x, y_i.i, "wtmean"))
  oos.ls <- map(responses, ~calc_ensemble(oos.ls, cv.ls, .x, y_i.i, "GLM_oos", ens.dir))
  
  saveRDS(fit.ls, glue("out/1_current/compiled/{y}_fit.rds"))
  saveRDS(oos.ls, glue("out/1_current/compiled/{y}_oos.rds"))
  
  
  
  # . null ------------------------------------------------------------------
  
  fit.ls <- readRDS(glue("out/1_current/compiled/{y}_fit.rds"))
  oos.ls <- readRDS(glue("out/1_current/compiled/{y}_oos.rds"))
  
  null.ls <- map(responses, ~calc_null(fit.ls, .x))
  fit.ls <- map(null.ls, ~.x$obs.df)
  oos.ls <- map2(oos.ls, null.ls,
                 ~left_join(.x |> mutate(yday=yday(date)), .y$yday.df) |> select(-yday)) |>
    map2(.x=_, fit.ls, ~bind_cols(.x, .y |> select(contains("nullGrand")) |> slice_head(n=1)))
  
  saveRDS(fit.ls, glue("out/1_current/compiled/{y}_fit.rds"))
  saveRDS(oos.ls, glue("out/1_current/compiled/{y}_oos.rds"))
  saveRDS(map(null.ls, ~.x$yday.df), glue("out/1_current/compiled/{y}_null.rds"))
  
}

closeAllConnections()







# setup -------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", "furrr")
lapply(pkgs, library, character.only=T)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
source("code/00_fn.R")

y_i <- bind_rows(read_csv("data/i_hab.csv") |> arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv") |> arrange(abbr) |> mutate(type="tox"))
y_resp <- c("Alsp", "PSP", "Disp", "DSP", "Pssp", "Psde", "Psse", "ASP", "Kami")

mod_i <- tibble(levels=c("nullGrand", "null4wk", "nullAuto", "perfect",
                         "ens", "ensLogitMn", "ensGLM", "ensGLM2",
                         "Ridge", "MARS", "NN", 
                         "RF", "Boost",
                         "HBL1"),
                labels=c("Null (int.)", "Null (\u00B12wk avg)", "Null (auto)", "perfect", 
                         "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "Ens-Ridge2", 
                         "Ridge", "MARS", "NN",
                         "RF", "XGB",
                         "HB"))
mod_cols <- c(rep("grey", 3), "grey30",
              rep("grey30", 4),
              "#b2df8a", "#33a02c", "#ff7f00", 
              "#cab2d6", "#6a3d9a",
              "#1f78b4") |>
  setNames(mod_i$labels)
d_i <- tibble(f=dir("out/1_current/model_fits")) |>
  mutate(Regional=grepl("Avg1", f),
         EnvInt=grepl("Xf1", f),
         AutoInt=grepl("XN1", f),
         EnvDelta=grepl("Del1", f),
         covSet=paste0("d", str_split_fixed(f, "-", 2)[,1])) |>
  arrange(Regional, EnvDelta, AutoInt, EnvInt) |>
  mutate(covSet_reorder=factor(covSet, levels=unique(covSet)))

fit.ls <- dirf("out/1_current/compiled", "_fit.rds") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)
oos.ls <- dirf("out/1_current/compiled", "_oos.rds") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)
fit.ls$alert_L <- fit.ls$alert |> 
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3),
         nullAuto_A1=if_else(prevAlert=="A0", 1e-3, 1-1e-3)) |>
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") |> 
  mutate(model=str_split_fixed(run, "_", 3)[,1], 
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model), 
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) |> 
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(paste0("d", 1:16), "ens", "ensLogitMn", "ensGLM", "ensGLM2", 
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  left_join(d_i |> select(-f)) |>
  mutate(covSet=factor(covSet, levels=c(levels(d_i$covSet_reorder), 
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2", 
                                        "null4wk", "nullAuto", "nullGrand", "perfect")),
         y=factor(y, levels=y_resp))

oos.ls$alert_L <- oos.ls$alert |> 
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3),
         nullAuto_A1=if_else(prevAlert=="A0", 1e-3, 1-1e-3)) |>
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") |> 
  mutate(model=str_split_fixed(run, "_", 3)[,1], 
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model), 
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) |> 
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(paste0("d", 1:16), "ens", "ensLogitMn", "ensGLM", "ensGLM2", 
                                        "null4wk", "nullAuto", "nullGrand", "perfect"))) |>
  left_join(d_i |> select(-f)) |>
  mutate(covSet=factor(covSet, levels=c(levels(d_i$covSet_reorder), 
                                        "ens", "ensLogitMn", "ensGLM", "ensGLM2", 
                                        "null4wk", "nullAuto", "nullGrand", "perfect")),
         y=factor(y, levels=y_resp))





# Threshold analysis ------------------------------------------------------

gc()
thresh.fit <- list(
  compute_thresholds(fit.ls$alert_L |> filter(!grepl("Null", model)), 
                     0.01, 0.9, 0.03, byPrevAlert=T, cores=20),
  compute_thresholds(fit.ls$alert_L |> filter(grepl("Null", model)), 
                     0.01, 0.9, 0.03, byPrevAlert=F, cores=20)
)
gc()
opt.F1 <- list(
  thresh.fit[[1]] |> filter(!is.na(F1)) |>
    group_by(y, model, PCA, covSet, prevAlert) |>
    arrange(desc(F1)) |> slice_head(n=1) |> ungroup() |>
    select(y, model, PCA, covSet, thresh, F1, precision, recall, accuracy, prevAlert) |>
    rename(optF1=thresh, F1_fit=F1, F1_precision=precision, F1_recall=recall, F1_accuracy=accuracy),
  thresh.fit[[2]] |> filter(!is.na(F1)) |>
    group_by(y, model, PCA, covSet) |>
    arrange(desc(F1)) |> slice_head(n=1) |> ungroup() |>
    select(y, model, PCA, covSet, thresh, F1, precision, recall, accuracy) |>
    rename(optF1=thresh, F1_fit=F1, F1_precision=precision, F1_recall=recall, F1_accuracy=accuracy)
)
opt.F2 <- list(
  thresh.fit[[1]] |> filter(!is.na(F2)) |>
    group_by(y, model, PCA, covSet, prevAlert) |>
    arrange(desc(F2)) |> slice_head(n=1) |> ungroup() |>
    select(y, model, PCA, covSet, thresh, F2, precision, recall, accuracy, prevAlert) |>
    rename(optF2=thresh, F2_fit=F2, F2_precision=precision, F2_recall=recall, F2_accuracy=accuracy),
  thresh.fit[[2]] |> filter(!is.na(F2)) |>
    group_by(y, model, PCA, covSet) |>
    arrange(desc(F2)) |> slice_head(n=1) |> ungroup() |>
    select(y, model, PCA, covSet, thresh, F2, precision, recall, accuracy) |>
    rename(optF2=thresh, F2_fit=F2, F2_precision=precision, F2_recall=recall, F2_accuracy=accuracy)
)
oos.ls$alert_L <- bind_rows(
  oos.ls$alert_L |>
    filter(!grepl("Null", model)) |>
    left_join(opt.F1[[1]] |> select(-F1_fit)) |>
    left_join(opt.F2[[1]] |> select(-F2_fit)),
  oos.ls$alert_L |>
    filter(grepl("Null", model)) |>
    left_join(opt.F1[[2]] |> select(-F1_fit)) |>
    left_join(opt.F2[[2]] |> select(-F2_fit))
) |>
  mutate(predF1=factor(if_else(prA1 > optF1, "A1", "A0"), levels=c("A0", "A1")),
         predF2=factor(if_else(prA1 > optF2, "A1", "A0"), levels=c("A0", "A1")))


saveRDS(fit.ls$alert_L, "out/1_current/clean/out_fit.rds")
saveRDS(oos.ls$alert_L, "out/1_current/clean/out_oos.rds")



library(kerneval)
rank.df <-  oos.ls$alert_L |> 
  select(y, model, covSet, PCA, alert, prA1) |>
  filter(!grepl("perfect|auto", model)) |>
  na.omit() |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y) |>
  mutate(rank=min_rank(desc(AUCNPR)),
         .metric="PR-AUC") |>
  rename(.estimate=AUCNPR) |> select(-AUCPR) |>
  bind_rows(oos.ls$alert_L |> 
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              roc_auc(prA1, truth=alert, event_level="second") |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="ROC-AUC") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              select(y, covSet, PCA, model, obsid, alert, prA1) |>
              filter(!is.na(prA1)) |>
              pivot_wider(names_from="alert", values_from="prA1") |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=schoenr(density(A0, na.rm=T), density(A1, na.rm=T))) |>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Schoener's D")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              f_meas(predF1, truth=alert, beta=1, event_level="second") |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="F1") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              accuracy(predF1, truth=alert, event_level="second") |>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="Accuracy (F1)") |>
              select(-.estimator)) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(predF1=="A1"))|>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Precision: TP/(TP+FP) (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(alert=="A1"))|>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="Recall: TP/(TP+FN) (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A0")/n())|>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="FPR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=sum(predF1=="A1" & alert=="A1")/n())|>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="TPR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=sum(predF1=="A0" & alert=="A1")/n())|>
              group_by(y) |>
              mutate(rank=min_rank(.estimate),
                     .metric="FNR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              group_by(y, model, PCA, covSet) |>
              summarise(.estimate=sum(predF1=="A0" & alert=="A0")/n())|>
              group_by(y) |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="TNR (F1)")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="vz", y) |>
              rename(.estimate=R2) |>
              mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
              na.omit() |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-VZ")) |>
  bind_rows(oos.ls$alert_L |>
              filter(!grepl("perfect|auto", model)) |>
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1),
                     alert=as.numeric(alert=="A1")) |>
              calc_R2(type="mf", y) |>
              rename(.estimate=R2) |>
              na.omit() |>
              mutate(rank=min_rank(desc(.estimate)),
                     .metric="R2-MF"))

saveRDS(rank.df, "out/1_current/clean/rank_oos.rds")


