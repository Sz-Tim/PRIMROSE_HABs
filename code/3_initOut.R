# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model analysis



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
                         "ensHB", "ensRF", "ensRF2", 
                         "Ridge", "MARS", "NN", 
                         "RF", "Boost",
                         "HBL1"),
                labels=c("Null (int.)", "Null (\u00B12wk avg)", "Null (auto)", "perfect", 
                         "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "Ens-Ridge2", 
                         "Ens-HB", "Ens-RF", "Ens-RF2", 
                         "Ridge", "MARS", "NN",
                         "RF", "XGB",
                         "HB"))
# labels=c("Null (intercept)", "Null (\u00B12wk avg)", "Null (auto)", "perfect", 
#          "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "Ens-Ridge2", 
#          "Ens-HB", "Ens-RF", "Ens-RF2", 
#          "Ridge GLM", "MARS", "NeuralNetwork",
#          "RandForest", "XGBoost",
#          "Bayes"))
mod_cols <- c(rep("grey", 3), "grey30",
              rep("grey30", 7),
              "#b2df8a", "#33a02c", "#ff7f00", 
              "#cab2d6", "#6a3d9a",
              "#1f78b4") |>
  setNames(mod_i$labels)
d_i <- tibble(f=dir("out/model_fits")) |>
  mutate(Regional=grepl("Avg1", f),
         EnvInt=grepl("Xf1", f),
         AutoInt=grepl("XN1", f),
         EnvDelta=grepl("Del1", f),
         covSet=paste0("d", str_split_fixed(f, "-", 2)[,1])) |>
  arrange(Regional, EnvDelta, AutoInt, EnvInt) |>
  mutate(covSet_reorder=factor(covSet, levels=unique(covSet)))

fit.ls <- dirf("out/compiled", "_fit.rds") |>
  map(~readRDS(.x)) |> list_transpose() |> map(bind_rows)
oos.ls <- dirf("out/compiled", "_oos.rds") |>
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
                     0.01, 0.9, 0.05, byPrevAlert=T, cores=20),
  compute_thresholds(fit.ls$alert_L |> filter(grepl("Null", model)), 
                     0.01, 0.9, 0.05, byPrevAlert=F, cores=20)
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


saveRDS(fit.ls$alert_L, "out/clean/out_fit.rds")
saveRDS(oos.ls$alert_L, "out/clean/out_oos.rds")



# metrics -----------------------------------------------------------------

# Metric summaries
oos.ls$alert_L |> 
  filter(model != "perfect") |>
  na.omit() |>
  group_by(y, model, covSet, PCA) |>
  roc_auc(prA1, truth=alert, event_level="second") |>
  group_by(y) |> arrange(desc(.estimate)) |> slice_head(n=1) |>
  ungroup() |> arrange(desc(.estimate))
oos.ls$alert_L |> 
  select(y, model, covSet, PCA, alert, prA1) |>
  filter(model != "perfect") |>
  na.omit() |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y) |> arrange(model) |> 
  mutate(delta=AUCNPR-first(AUCNPR), deltaPct=100*delta/(1-first(AUCNPR))) |>
  group_by(y) |> arrange(desc(AUCNPR)) |> slice_head(n=1) |>
  ungroup() |> arrange(desc(AUCNPR))
oos.ls$alert_L |> 
  select(y, model, covSet, PCA, alert, prA1, prevAlert) |>
  filter(model != "perfect") |>
  filter(model != "Null (auto)") |>
  na.omit() |>
  find_AUCPR_min(y, prevAlert) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y, prevAlert) |> arrange(desc(AUCNPR)) |> slice_head(n=1) |>
  ungroup() |> arrange(desc(AUCNPR))

oos.ls$alert_L |> 
  filter(grepl("Null|Ens", model)) |>
  select(y, model, covSet, alert, prA1) |>
  na.omit() |>
  find_AUCPR_min(y) |>
  filter(grepl("wk|Ens", model)) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  group_by(y) |> arrange(model) |> 
  mutate(delta=AUCNPR-first(AUCNPR), deltaPct=100*delta/(1-first(AUCNPR))) |>
  ungroup() |> filter(!grepl("Null", model)) |>
  arrange(y, desc(AUCNPR)) |>
  print(n=50)


# Plots
oos.ls$alert_L |> 
  filter(grepl("^d|null4wk", covSet)) |>
  mutate(covSet=as.character(covSet)) |>
  group_by(y, model, covSet, PCA) |>
  roc_auc(prA1, truth=alert, event_level="second") |>
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0.5, 1) + labs(x="", y="ROC-AUC (oos)") + 
  facet_grid(.~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L |> 
  select(y, model, covSet, PCA, alert, prA1) |>
  filter(model != "perfect") |>
  filter(model != "Null (auto)") |>
  na.omit() |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  filter(model != "Null (intercept)") |>
  ggplot(aes(covSet, AUCNPR, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Normalized PR-AUC (oos)") + 
  facet_grid(.~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L |> 
  select(y, model, covSet, PCA, alert, prA1, prevAlert) |>
  filter(model != "perfect") |>
  filter(model != "Null (auto)") |>
  na.omit() |>
  find_AUCPR_min(y, prevAlert) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  filter(model != "Null (intercept)") |>
  ggplot(aes(covSet, AUCNPR, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Normalized PR-AUC (oos)") + 
  facet_grid(prevAlert~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))


oos.ls$alert_L |> 
  filter(!is.na(prA1)) |>
  group_by(y, model, PCA, covSet) |>
  summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(predF1=="A1")) |>
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Precision: Pr(A1 | predA1) (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L |> 
  filter(!is.na(prA1)) |>
  group_by(y, model, PCA, covSet) |>
  summarise(.estimate=sum(predF1=="A1" & alert=="A1")/sum(alert=="A1")) |>
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Recall: Pr(predA1 | A1) (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L |> 
  group_by(y, model, PCA, covSet) |>
  f_meas(predF1, truth=alert, beta=1, event_level="second") |>
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="F1 (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L |> 
  group_by(y, model, PCA, covSet) |>
  f_meas(predF2, truth=alert, beta=2, event_level="second") |>
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="F2 (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))



oos.ls$alert_L |>
  filter(!grepl("perfect", model)) |>
  filter(grepl("intercept|2wk|Ens", model)) |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1),
         alert=as.numeric(alert=="A1")) |>
  group_by(y, model, PCA, covSet) |>
  calc_R2(type="vz") |>
  filter(!is.na(R2)) |>
  ggplot(aes(covSet, R2, colour=model, linetype=PCA, shape=PCA)) + 
  geom_point() +
  geom_line(aes(group=paste(model, PCA)), linewidth=0.8) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  labs(x="Model", y="VZ R2 (oos)") + 
  facet_grid(.~y) + 
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")

oos.ls$alert_L |>
  filter(!grepl("perfect", model)) |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1),
         alert=as.numeric(alert=="A1")) |>
  calc_R2(type="mf", y) |>
  filter(model != "Null (intercept)") |>
  filter(!is.na(R2)) |>
  ggplot(aes(covSet, R2, colour=model, linetype=PCA, shape=PCA)) + 
  geom_point() +
  geom_line(aes(group=paste(model, PCA)), linewidth=0.8) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, NA) + labs(x="Model", y="McFaddens R2 (oos)") + 
  facet_grid(.~y) +
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")
oos.ls$alert_L |>
  filter(!grepl("perfect", model)) |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1),
         alert=as.numeric(alert=="A1")) |>
  calc_R2(type="mf", y, siteid) |>
  filter(model != "Null (intercept)") |>
  group_by(y, model, PCA, covSet) |>
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n())) |>
  filter(!is.na(R2_mn)) |>
  ggplot(aes(covSet, R2_mn, colour=model)) + 
  geom_point(aes(shape=PCA)) + 
  geom_line(aes(linetype=PCA, group=paste(model, PCA)), linewidth=0.8) + 
  # geom_errorbar(aes(ymin=R2_mn-2*se, ymax=R2_mn+2*se), width=0.5) +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  labs(x="Model", y="McFaddens R2 mn, 2se among sites (oos)") + 
  facet_grid(.~y) + ylim(0, NA) + 
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")
oos.ls$alert_L |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) |>
  group_by(y, siteid) |>
  arrange(y, model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  ggplot(aes(covSet, R2, fill=model)) + 
  geom_boxplot() +
  scale_fill_manual(values=mod_cols) +
  labs(x="Model", y="McFaddens R2 mn, 2se among sites (oos)") + 
  ylim(-1, 1) + 
  facet_wrap(~y) +
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")




oos_precCovSet <- oos.ls$alert_L |> 
  filter(grepl("^d", covSet)) |>
  na.omit() |>
  group_by(y, model, PCA, Regional, EnvInt, AutoInt, EnvDelta) |>
  average_precision(prA1, truth=alert, event_level="second")

oos_precCovSet |>
  group_by(y, model, Regional, EnvInt, AutoInt, EnvDelta) |>
  arrange(PCA) |>
  summarise(deltaVsF=last(.estimate)-first(.estimate)) |>
  ungroup() |>
  ggplot(aes(deltaVsF, model, colour=model)) + 
  geom_vline(xintercept=0) + 
  geom_boxplot() + 
  scale_colour_manual(values=mod_cols, guide="none") + 
  facet_wrap(~y) + xlim(-0.5, 0.5) +
  xlab("AUCPR: PCA - noPCA")
oos_precCovSet |>
  group_by(y, model, PCA, EnvInt, AutoInt, EnvDelta) |>
  arrange(Regional) |>
  summarise(deltaVsF=last(.estimate)-first(.estimate)) |>
  ungroup() |>
  ggplot(aes(deltaVsF, model, colour=model)) + 
  geom_vline(xintercept=0) + 
  geom_boxplot() + 
  scale_colour_manual(values=mod_cols, guide="none") + 
  facet_wrap(~y) + xlim(-0.5, 0.5) +
  xlab("AUCPR: Regional - noRegional")
oos_precCovSet |>
  group_by(y, model, PCA, Regional, AutoInt, EnvDelta) |>
  arrange(EnvInt) |>
  summarise(deltaVsF=last(.estimate)-first(.estimate)) |>
  ungroup() |>
  ggplot(aes(deltaVsF, model, colour=model)) + 
  geom_vline(xintercept=0) + 
  geom_boxplot() + 
  scale_colour_manual(values=mod_cols, guide="none") + 
  facet_wrap(~y) + xlim(-0.5, 0.5) +
  xlab("AUCPR: EnvInt - noEnvInt")
oos_precCovSet |>
  group_by(y, model, PCA, Regional, EnvInt, EnvDelta) |>
  arrange(AutoInt) |>
  summarise(deltaVsF=last(.estimate)-first(.estimate)) |>
  ungroup() |>
  ggplot(aes(deltaVsF, model, colour=model)) + 
  geom_vline(xintercept=0) + 
  geom_boxplot() + 
  scale_colour_manual(values=mod_cols, guide="none") + 
  facet_wrap(~y) + xlim(-0.5, 0.5) +
  xlab("AUCPR: AutoInt - noAutoInt")
oos_precCovSet |>
  group_by(y, model, PCA, Regional, EnvInt, AutoInt) |>
  arrange(EnvDelta) |>
  summarise(deltaVsF=last(.estimate)-first(.estimate)) |>
  ungroup() |>
  ggplot(aes(deltaVsF, model, colour=model)) + 
  geom_vline(xintercept=0) + 
  geom_boxplot() + 
  scale_colour_manual(values=mod_cols, guide="none") + 
  facet_wrap(~y) + xlim(-0.5, 0.5) +
  xlab("AUCPR: EnvDelta - noEnvDelta")


oos.ls$alert_L |>
  filter(grepl("d", covSet)) |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  group_by(y, model, PCA, covSet) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) |>
  group_by(y) |>
  arrange(y, model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  filter(!is.na(R2)) |>
  mutate(R2=if_else(R2 < 0, 0, R2)) |>
  ggplot(aes(covSet, R2, colour=PCA, group=PCA)) + geom_point() + geom_line() +
  facet_grid(y~model)





oos.ls$alert_L |> 
  filter(grepl("^d", covSet)) |>
  na.omit() |>
  group_by(y, model, covSet, PCA) |>
  average_precision(prA1, truth=alert, event_level="second") |>
  ggplot(aes(covSet, .estimate, colour=model, group=model)) + 
  geom_point() + geom_line() +
  scale_colour_manual(values=mod_cols) +
  facet_grid(PCA~y) + labs(x="", y="PR-AUC")
oos.ls$alert_L |> 
  filter(grepl("^d", covSet)) |>
  na.omit() |>
  group_by(y, model, covSet, PCA) |>
  average_precision(prA1, truth=alert, event_level="second") |>
  ggplot(aes(covSet, .estimate, colour=PCA, group=PCA)) + 
  geom_point(size=2, alpha=0.5) + geom_line() +
  facet_grid(model~y) + labs(x="", y="PR-AUC")

oos.ls$alert_L |> 
  filter(grepl("^d", covSet)) |>
  na.omit() |>
  group_by(y, model, PCA, covSet) |>
  f_meas(predF, truth=alert, event_level="second") |>
  ggplot(aes(covSet, .estimate, colour=model, group=model)) + 
  geom_point() + geom_line() +
  scale_colour_manual(values=mod_cols) +
  facet_grid(PCA~y) + labs(x="", y="F1")


oos.ls$alert_L |>
  filter(model != "perfect") |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  group_by(y, model, PCA, covSet) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) |>
  group_by(y) |>
  arrange(model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  arrange(desc(R2)) |>
  slice_head(n=1) |>
  ungroup() |>
  arrange(desc(R2))

oos.ls$alert_L |>
  filter(grepl("intercept|wk|Ens", model)) |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  group_by(y, model, PCA, covSet) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) |>
  group_by(y) |>
  arrange(model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  group_by(y) |>
  arrange(model) |>
  mutate(delta=R2-first(R2)) |>
  # filter(model=="Ens") |>
  ungroup() |>
  arrange(desc(delta))


oos.ls$alert_L |>
  filter(model != "perfect") |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  # filter(PCA != "lnN") |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) |>
  group_by(y, siteid) |>
  arrange(y, model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  group_by(y, model, PCA, covSet) |>
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n()), nA1=sum(nA1), mnPA1=mean(pA1)) |>
  ungroup() |>
  filter(!is.na(R2_mn)) |>
  group_by(y) |>
  arrange(R2_mn) |> 
  slice_tail(n=1) |>
  arrange(desc(R2_mn))

oos.ls$alert_L |>
  filter(grepl("intercept|wk|Ens", model)) |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  # filter(PCA != "lnN") |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) |>
  group_by(y, siteid) |>
  arrange(y, model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  group_by(y, model, PCA, covSet) |>
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n()), nA1=sum(nA1), mnPA1=mean(pA1)) |>
  filter(model != "Null (intercept)") |>
  group_by(y) |>
  arrange(model) |>
  mutate(delta=R2_mn-first(R2_mn)) |>
  # filter(model=="Ens-Ridge") |>
  ungroup() |>
  arrange(desc(delta))


oos.ls$alert_L |>
  filter(model != "perfect") |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  # filter(PCA != "lnN") |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) |>
  group_by(y, siteid) |>
  arrange(y, model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  group_by(y, model, PCA, covSet) |>
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n()), nA1=sum(nA1), mnPA1=mean(pA1)) |>
  ungroup() |>
  filter(!is.na(R2_mn)) |>
  ggplot(aes(mnPA1, R2_mn)) + 
  geom_point() + facet_wrap(~model)


site.df <- bind_rows(
  readRDS("data/site_hab_df.rds") |> mutate(type="hab"),
  readRDS("data/site_tox_df.rds") |> mutate(type="tox")
) |>
  select(type, siteid, sin, lon, lat, fetch, openBearing)
library(sf)
coast <- st_read("data/northAtlantic_footprint.gpkg") |>
  st_crop(st_bbox(c(xmin=min(site.df$lon)-1e4, xmax=max(site.df$lon)+1e4,
                  ymin=min(site.df$lat)-1e4, ymax=max(site.df$lat)+1e4), crs=27700))
siteR2 <- oos.ls$alert_L |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) |>
  group_by(y, siteid) |>
  arrange(y, model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  mutate(type=y_i$type[match(y, y_i$abbr)]) |>
  left_join(site.df)

siteR2 |>
  filter(pA1 > 0 & pA1 < 1) |>
  ggplot(aes(pA1, R2, colour=model)) + ylim(-1,1) + 
  scale_colour_manual(values=mod_cols) +
  geom_point(alpha=0.6, shape=1) + 
  stat_smooth(se=F, method="lm") + 
  facet_wrap(~y) +
  scale_x_continuous(trans="logit") +
  theme(legend.position="none")
siteR2 |> 
  filter(covSet=="ensGLM") |>
  ggplot(aes(pA1, R2)) + 
  geom_point() +
  facet_wrap(~y)
siteR2 |>
  filter(covSet=="ensGLM") |>
  ggplot(aes(fetch, R2, colour=pA1, group=siteid)) +
  geom_point() + 
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
siteR2 |>
  filter(covSet=="ensGLM") |>
  ggplot(aes(openBearing, R2)) +
  geom_point() + stat_smooth(method="gam", formula=y~s(x,bs="cc")) +
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
siteR2 |>
  filter(covSet=="ensGLM") |>
  ggplot() + 
  geom_sf(data=coast, colour=NA) + 
  geom_point(aes(lon, lat, colour=pA1), shape=1, size=2) + 
  scale_colour_viridis_c() +
  facet_wrap(~y, nrow=2) + theme_classic()
siteR2 |>
  filter(covSet=="ensGLM") |>
  ggplot() + 
  geom_sf(data=coast, colour=NA) + 
  geom_point(aes(lon, lat, colour=R2), shape=1, size=2) + 
  scale_colour_viridis_c(option="turbo", limits=c(0,1)) +
  facet_wrap(~y, nrow=2) + theme_classic()

siteAUC <- oos.ls$alert_L |>
  filter(!grepl("intercept", model)) |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(pA1=mean(alert=="A1")) |>
  ungroup() |>
  filter(pA1 > 0) |>
  group_by(y, model, PCA, covSet, siteid, pA1) |>
  roc_auc(prA1, truth=alert, event_level="second") |>
  mutate(type=y_i$type[match(y, y_i$abbr)]) |>
  left_join(site.df)
siteAUC |>
  ggplot(aes(covSet, .estimate, fill=model)) + ylim(0,1) + 
  scale_fill_manual(values=mod_cols) +
  geom_boxplot() + 
  facet_wrap(~y) +
  theme(legend.position="none")
siteAUC |>
  ggplot(aes(pA1, .estimate, colour=model)) + ylim(0,1) + 
  scale_colour_manual(values=mod_cols) +
  geom_point(alpha=0.6, shape=1) + 
  stat_smooth(se=F, method="lm") + 
  facet_wrap(~y) +
  scale_x_continuous(trans="logit") +
  theme(legend.position="none")
siteAUC |>
  filter(covSet=="ensGLM") |>
  ggplot(aes(fetch, .estimate)) +
  geom_point(aes(colour=pA1)) + stat_smooth(method="loess", se=F, span=2) +
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
siteAUC |>
  filter(covSet=="ensGLM") |>
  ggplot(aes(openBearing, .estimate)) +
  geom_point() + stat_smooth(method="gam", formula=y~s(x,bs="cc")) +
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
siteAUC |>
  filter(covSet=="ensGLM") |>
  ggplot() + 
  geom_sf(data=coast, colour=NA) + 
  geom_jitter(aes(lon, lat, colour=.estimate), size=2, width=1e4, height=1e4) + 
  scale_colour_viridis_c("ROC-AUC") +
  facet_wrap(~y, nrow=2) + theme_classic()
siteAUC |> 
  group_by(y, siteid) |>
  summarise(mnAUC=median(.estimate, na.rm=T), pA1=first(pA1)) |>
  ggplot(aes(pA1, mnAUC)) + 
  geom_point() +
  facet_wrap(~y)

sitePRAUC <- oos.ls$alert_L |> 
  select(y, model, covSet, PCA, siteid, alert, prA1) |>
  filter(model != "perfect") |>
  na.omit() |>
  find_AUCPR_min(y, siteid) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         .estimate=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  mutate(type=y_i$type[match(y, y_i$abbr)]) |>
  left_join(site.df)
sitePRAUC <- oos.ls$alert_L |>
  filter(!is.na(prA1)) |>
  find_AUCPR_min(y, siteid) |>
  filter(!grepl("intercept", model)) |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(pA1=mean(alert=="A1")) |>
  ungroup() |>
  filter(pA1 > 0) |>
  group_by(y, model, PCA, covSet, siteid, pA1, AUCPR_min) |>
  average_precision(prA1, truth=alert, event_level="second") |>
  mutate(.estimate=(.estimate-AUCPR_min)/(1-AUCPR_min),
         type=y_i$type[match(y, y_i$abbr)]) |>
  left_join(site.df)
sitePRAUC |>
  ggplot(aes(pA1, .estimate, colour=model)) + ylim(0,1) + 
  scale_colour_manual(values=mod_cols) +
  geom_point(alpha=0.6, shape=1) + 
  stat_smooth(se=F, method="lm") + 
  facet_wrap(~y) +
  scale_x_continuous(trans="logit") +
  theme(legend.position="none")
sitePRAUC |>
  filter(covSet=="ensGLM") |>
  ggplot(aes(fetch, .estimate)) +
  geom_point() + 
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
sitePRAUC |>
  filter(covSet=="ensGLM") |>
  ggplot(aes(openBearing, .estimate)) +
  geom_point() + stat_smooth(method="gam", formula=y~s(x,bs="cc")) +
  facet_wrap(~y, scales="free_x")
sitePRAUC |>
  filter(covSet=="ensGLM") |>
  ggplot() + 
  geom_sf(data=coast, colour=NA) + 
  geom_point(aes(lon, lat, colour=.estimate), shape=1, size=2) + 
  geom_point(aes(lon, lat, colour=.estimate), shape=1, size=1.5) + 
  geom_point(aes(lon, lat, colour=.estimate), shape=1, size=1) + 
  scale_colour_viridis_c(option="turbo", limits=c(0,1)) +
  facet_wrap(~y, nrow=2) + theme_classic()

sitePRAUC |>
  ggplot(aes(siteid, .estimate, colour=pA1, group=siteid)) +
  geom_boxplot() + 
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
sitePRAUC |> 
  group_by(y, siteid) |>
  summarise(mnAUC=median(.estimate, na.rm=T), pA1=first(pA1)) |>
  ggplot(aes(pA1, mnAUC)) + 
  geom_point() +
  facet_wrap(~y)



oos.ls$alert_L |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) |>
  # filter(PCA != "lnN") |>
  group_by(y, model, PCA, covSet, siteid) |>
  mutate(alert=as.numeric(alert=="A1")) |>
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) |>
  group_by(y, siteid) |>
  arrange(y, model) |>
  mutate(R2=1 - LL/first(LL)) |>
  filter(model != "Null (intercept)") |>
  mutate(anyA1=pA1 > 0) |>
  ggplot(aes(anyA1, R2, fill=model)) + ylim(-1,1) + 
  scale_fill_manual(values=mod_cols) +
  geom_boxplot() +
  facet_wrap(~y) +
  theme(legend.position="none")



oos.ls$alert_L |> 
  filter(grepl("2wk|Ens-R", model)) |>
  # na.omit() |>
  group_by(y, model, PCA, covSet) |> 
  pr_curve(prA1, truth=alert, event_level="second") |> 
  filter(recall > 0) |>
  ggplot(aes(recall, precision, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y) + ggtitle("oos")
oos.ls$alert_L |> 
  filter(grepl("2wk|intercept|Ens|perfect", model)) |>
  group_by(y, model, PCA, covSet) |> 
  gain_curve(prA1, truth=alert, event_level="second") |> 
  ggplot(aes(.percent_tested, .percent_found, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)
oos.ls$alert_L |> 
  filter(grepl("2wk|Ens-R", model)) |>
  # na.omit() |>
  group_by(y, model, PCA, covSet) |> 
  roc_curve(prA1, truth=alert, event_level="second") |> 
  ggplot(aes(1-specificity, sensitivity, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + geom_abline() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)

fit.ls$alert_L |> 
  # na.omit() |>
  filter(grepl("2wk|Ens", model)) |>
  group_by(y, model, PCA, covSet) |> 
  pr_curve(prA1, truth=alert, event_level="second") |> 
  filter(recall > 0) |>
  ggplot(aes(recall, precision, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y) + ggtitle("fit")
fit.ls$alert_L |> 
  filter(grepl("2wk|intercept|Ens|perfect", model)) |>
  group_by(y, model, PCA, covSet) |> 
  gain_curve(prA1, truth=alert, event_level="second") |> 
  ggplot(aes(.percent_tested, .percent_found, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)
fit.ls$alert_L |> 
  filter(grepl("2wk|Ens", model)) |>
  # na.omit() |>
  group_by(y, model, PCA, covSet) |> 
  roc_curve(prA1, truth=alert, event_level="second") |> 
  ggplot(aes(1-specificity, sensitivity, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + geom_abline() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)



# ranks -------------------------------------------------------------------


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
              select(y, covSet, PCA, model, obsid, alert, prA1) %>%
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

saveRDS(rank.df, "out/clean/rank_oos.rds")

rank.df |>
  filter(.metric!="roc_auc") |>
  group_by(model) |>
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) |>
  arrange(meanRank)
rank.df |> group_by(model, .metric) |>
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) |>
  arrange(.metric, meanRank) |>
  filter(grepl("Ens", model)) |>
  print(n=21)

rank.df |> group_by(covSet) |>
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) |>
  arrange(meanRank)
rank.df |> group_by(covSet, PCA) |>
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) |>
  arrange(meanRank)


ggplot(rank.df, aes(.metric, rank, fill=model)) + 
  geom_boxplot() +
  scale_fill_manual(values=mod_cols) + 
  labs(x="", y="Rank") +
  facet_wrap(~.metric, scales="free_x", nrow=1) + 
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom")

rank.df |> 
  filter(.metric != "roc_auc") |>
  filter(grepl("^d", covSet)) |>
  ggplot(aes(.metric, rank, fill=covSet)) + 
  geom_boxplot() +
  scale_fill_viridis_d() + 
  labs(x="", y="Rank") + facet_wrap(~model)

rank.df |>
  filter(grepl("AUC|R2", .metric)) |>
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA, group=paste(model, PCA))) + 
  geom_point() + geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + facet_grid(.metric~y) 
  # ylim(0, NA) + facet_grid(.metric~y, scales="free_y")
rank.df |>
  filter(grepl("4wk|ens", covSet)) |>
  ggplot(aes(covSet, .estimate, colour=model)) + 
  geom_point() + 
  scale_colour_manual(values=mod_cols) + 
  scale_shape_manual(values=c(19, 1)) +
  facet_grid(.metric~y, scales="free_y")
rank.df |>
  filter(grepl("4wk|ensGLM|ensHB", covSet)) |>
  filter(grepl("AUC|R2|D|Accuracy", .metric)) |>
  # filter(grepl("Grand|4wk|ensG", covSet)) |>
  ggplot(aes(model, .estimate, colour=y, group=y)) + 
  geom_point() + geom_line() +
  scale_colour_brewer("", type="qual", palette=3) +
  ylab("Value") +
  facet_wrap(~.metric, nrow=1) + ylim(NA,1) +
  # facet_wrap(~.metric, nrow=1, scales="free_y") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        legend.position="bottom")

rank.df |>
  filter(grepl("Bayes", model)) |>
  mutate(covSet=factor(covSet, levels=paste0("d", 1:16))) |>
  ggplot(aes(covSet, .estimate, colour=y, shape=PCA, group=paste(model, PCA, y))) + 
  geom_point() + geom_line() +
  scale_shape_manual(values=c(19, 1)) +
  scale_colour_brewer("", type="qual", palette=3) +
  ylab("Value") +
  facet_wrap(~.metric, nrow=1, scales="free_y") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        legend.position="bottom")




# probabilities -----------------------------------------------------------

oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge", model)) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_boxplot(outlier.size=0.5, outlier.shape=1, outlier.alpha=0.5) +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_wrap(~y) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge", model)) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_boxplot(outlier.size=0.5, outlier.shape=1, outlier.alpha=0.5) +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_grid(prevAlert~y) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L |>
  filter(grepl("2wk|Ens", model)) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_violin(scale="width") +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_wrap(~y) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge", model)) |>
  filter(!grepl("intercept", model)) |>
  group_by(y, alert, model) |>
  get_intervals(prA1, "qi") |>
  ggplot(aes(med, xmin=L05, xmax=L95, model, colour=alert)) + 
  geom_linerange(aes(xmin=L25, xmax=L75), position=position_dodge(width=0.5), linewidth=2) + 
  geom_linerange(aes(xmin=L10, xmax=L90), position=position_dodge(width=0.5), linewidth=1) + 
  geom_pointrange(position=position_dodge(width=0.5)) + 
  scale_colour_manual("Truth", values=c("grey", "red3")) +
  facet_wrap(~y) +
  xlim(0,1)

oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge", model)) |>
  filter(!grepl("intercept", model)) |>
  group_by(y, prevAlert, alert, model) |>
  get_intervals(prA1, "qi") |>
  ggplot(aes(med, xmin=L05, xmax=L95, model, colour=alert)) + 
  geom_linerange(aes(xmin=L25, xmax=L75), position=position_dodge(width=0.5), linewidth=2) + 
  geom_linerange(aes(xmin=L10, xmax=L90), position=position_dodge(width=0.5), linewidth=1) + 
  geom_linerange(position=position_dodge(width=0.5), linewidth=0.25) + 
  geom_point(position=position_dodge(width=0.5), size=3, shape=3) + 
  scale_colour_manual("Truth", values=c("grey", "red3")) +
  facet_grid(prevAlert~y) +
  xlim(0,1)

oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge", model)) |>
  ggplot(aes(prA1, model, fill=alert)) +
  ggridges::geom_density_ridges(rel_min_height=0.01, scale=1, alpha=0.5) +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_grid(.~y) +
  xlim(0,1)
oos.ls$alert_L |>
  filter(grepl("2wk|Ens", model)) |>
  ggplot(aes(prA1, model, fill=alert)) +
  ggridges::geom_density_ridges(rel_min_height=0.01, scale=1, alpha=0.5) +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_grid(prevAlert~y) +
  xlim(0,1)


oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge2", model)) |>
  mutate(pr_round=round(prA1, 1)) |>
  group_by(y, model, pr_round) |>
  summarise(true_alertPr=mean(alert=="A1"), N=n()) |>
  ggplot(aes(pr_round, true_alertPr)) + 
  geom_point(aes(colour=model, size=N), alpha=0.5) + 
  # geom_line(aes(colour=model)) + 
  stat_smooth(aes(colour=model), method="lm", se=F, linewidth=0.5) + 
  geom_abline(linetype=2) + 
  scale_colour_manual(values=mod_cols) +
  facet_wrap(~y) + 
  xlim(0,1) + ylim(0,1) + 
  labs(x="Forecast", y="Prorportion of true alerts")

oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge2", model)) |>
  mutate(pr_round=round(prA1, 1)) |>
  group_by(y, model, prevAlert, pr_round) |>
  summarise(true_alertPr=mean(alert=="A1"), N=n()) |>
  ggplot(aes(pr_round, true_alertPr)) + 
  geom_point(aes(colour=model, size=N), alpha=0.5) + 
  # geom_line(aes(colour=model)) + 
  stat_smooth(aes(colour=model), method="lm", se=F, linewidth=0.5) + 
  geom_abline(linetype=2) + 
  scale_colour_manual(values=mod_cols) +
  facet_grid(prevAlert~y) + 
  xlim(0,1) + ylim(0,1) + 
  labs(x="Forecast", y="Prorportion of true alerts")




# Mosaics
library(ggmosaic)
oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge2", model)) |>
  rename(Obs=alert, Pred=predF1, Previous=prevAlert) |>
  ggplot() + 
  geom_mosaic(aes(x=product(Pred, Obs),
                  fill=Pred)) + 
  scale_fill_manual(values=c("grey", "red3")) +
  facet_grid(y~model) +
  theme(panel.grid=element_blank())

oos.ls$alert_L |>
  filter(grepl("2wk|Ens-Ridge2", model)) |>
  mutate(Observed=factor(alert, labels=c("None", "Alert")),
         Predicted=factor(predF1, labels=c("None", "Alert")),
         Previous=factor(prevAlert, labels=c("t-1: None", "t-1: Alert"))) |>
  ggplot() + 
  geom_mosaic(aes(x=product(Predicted, Observed),
                  fill=Predicted)) + 
  scale_fill_manual(values=c("grey", "red3")) +
  facet_grid(y~model*Previous) +
  theme(panel.grid=element_blank())





# variable importance -----------------------------------------------------

col_metadata <- c("obsid", "y", "date", "year", "yday", "siteid", "lon", "lat")
col_resp <- c("lnN", "tl", "alert")
col_cmems <- readRDS("data/cmems_vars.rds")
col_wrf <- readRDS("data/wrf_vars.rds")

varTypes <- list(
  spacetime=c("ydayCos", "ydaySin", "ydaySinXydayCos",
              "latz", "lonz", "lonzXlatz"),
  main=c(
    "fetch",
    "lnNWt1", "lnNAvg1", "prAlertAvg1", "alert1A1",
    "lnNWt2", "lnNAvg2", "prAlertAvg2", "alert2A1",
    "lnNPrevYr", "lnNAvgPrevYr", "prAlertPrevYr", "prAlertAvgPrevYr",
    col_cmems, col_wrf
  ),
  UV_interact=c(
    paste("UWkXfetch", grep("Dir[EW]", col_cmems, value=T), sep="X"),
    paste("VWkXfetch", grep("Dir[NS]", col_cmems, value=T), sep="X"),
    paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X"),
    paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir[EW]", col_wrf, value=T), sep="X")
  ),
  hab=c(outer(filter(y_i, type=="hab")$abbr, c("lnNAvg", "prA"), "paste0"))
)
all_covs$interact <- paste("lnNWt1", c(all_covs$main[-2]), sep="X")

enet.f <- dirf("out/model_fits", "_Ridge.rds", recursive=T)
rf.f <- dirf("out/model_fits", "_RF.rds", recursive=T)
nn.f <- dirf("out/model_fits", "_NN.rds", recursive=T)
xgb.f <- dirf("out/model_fits", "_Boost.rds", recursive=T)

vi.df <- bind_rows(
  enet.f |> 
    map_dfr(~readRDS(.x) |> tidy() |> mutate(f=.x)) |> 
    filter(term!="(Intercept)") |>
    group_by(f) |> mutate(Importance=abs(estimate)/max(abs(estimate))) |> ungroup() |>
    rename(Variable=term) |> select(-estimate) |>
    mutate(model="Ridge"), 
  rf.f |>
    map_dfr(~readRDS(.x) |> extract_fit_engine() |> vip::vi(scale=T) |> mutate(f=.x)) |>
    group_by(f) |> mutate(Importance=Importance/max(Importance)) |> ungroup() |>
    mutate(model="RF"),
  # nn.f |>
  #   map_dfr(~readRDS(.x) |> extract_fit_engine() |> vip::vi(scale=T) |> mutate(f=.x)) |>
  #   group_by(f) |> mutate(Importance=abs(Importance)/max(abs(Importance))) |> ungroup() |>
  #   mutate(model="NN"),
  xgb.f |>
    map_dfr(~readRDS(.x) |> extract_fit_engine() |> vip::vi(scale=T) |> mutate(f=.x)) |>
    group_by(f) |> mutate(Importance=abs(Importance)/max(abs(Importance))) |> ungroup() |>
    mutate(model="XGB")
) |>
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         y=str_split_fixed(str_split_fixed(f, "/", 4)[,4], "_", 3)[,1]) 

vi.df |>
  group_by(y, Variable) |>
  summarise(Importance_mean=mean(Importance)) |>
  group_by(y) |>
  arrange(desc(Importance_mean)) |>
  slice_head(n = 10) |>
  ggplot(aes(Importance_mean, Variable)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~y, scales="free_y", nrow=1)

vi.df |>
  mutate(varType=)