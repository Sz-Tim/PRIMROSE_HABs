# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model analysis



# setup -------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", "furrr")
lapply(pkgs, library, character.only=T)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
source("code/00_fn.R")

y_i <- bind_rows(read_csv("data/i_hab.csv") %>% arrange(abbr) %>% mutate(type="hab"),
                 read_csv("data/i_tox.csv") %>% arrange(abbr) %>% mutate(type="tox"))

mod_i <- tibble(levels=c("nullGrand", "null4wk", "nullAuto", "perfect",
                         "ens", "ensLogitMn", "ensGLM",
                         "ensHB", "ensRF1", 
                         "Ridge", "MARS", "NN", 
                         "RF", "Boost",
                         "HBL1"),
                labels=c("Null (grand)", "Null (\u00B12wk avg)", "Null (auto)", "perfect", 
                         "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", 
                         "Ens-HB", "Ens-RF", 
                         "Ridge GLM", "MARS", "NeuralNetwork",
                         "RandForest", "XGBoost",
                         "Bayes"))
mod_cols <- c(rep("grey", 3), "grey30",
              rep("black", 5),
              "#b2df8a", "#33a02c", "#ff7f00", 
              "#cab2d6", "#6a3d9a",
              "#1f78b4") %>%
  setNames(mod_i$labels)
d_i <- tibble(f=dir("out/model_fits")) %>%
  mutate(Spatial=if_else(grepl("Avg1", f), "L+R", "L"),
         EnvInt=if_else(grepl("Xf1", f), "EnvInt", "No EnvInt"),
         AutoInt=if_else(grepl("XN1", f), "AutoInt", "No AutoInt"),
         EnvDelta=if_else(grepl("Del1", f), "EnvDelta", "No EnvDelta"),
         covSet=paste0("d", str_split_fixed(f, "-", 2)[,1]))

fit.ls <- dirf("out/compiled", "_fit.rds") %>%
  map(~readRDS(.x)) %>% list_transpose() %>% map(bind_rows)
oos.ls <- dirf("out/compiled", "_oos.rds") %>%
  map(~readRDS(.x)) %>% list_transpose() %>% map(bind_rows)
fit.ls$alert_L <- fit.ls$alert %>% 
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3),
         nullAuto_A1=if_else(prevAlert=="A0", 1e-3, 1-1e-3)) %>%
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") %>% 
  mutate(model=str_split_fixed(run, "_", 3)[,1], 
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model), 
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) %>% 
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(paste0("d", 1:16), "ens", "ensLogitMn",
                                        "null4wk", "nullAuto", "nullGrand", "perfect")))# %>%
  # left_join(d_i %>% select(-f))
oos.ls$alert_L <- oos.ls$alert %>% 
  mutate(perfect_A1=if_else(alert=="A0", 1e-3, 1-1e-3),
         nullAuto_A1=if_else(prevAlert=="A0", 1e-3, 1-1e-3)) %>%
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") %>% 
  mutate(model=str_split_fixed(run, "_", 3)[,1], 
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=if_else(grepl("^d", model), 
                       str_split_fixed(str_remove(model, "PCA."), "\\.", 2)[,2],
                       str_remove(model, "PCA."))) %>% 
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         covSet=factor(covSet, levels=c(paste0("d", 1:16), "ens", "ensLogitMn",
                                        "null4wk", "nullAuto", "nullGrand", "perfect")))# %>%
  # left_join(d_i %>% select(-f))
  




# Threshold analysis ------------------------------------------------------

gc()
thresh.fit <- list(
  compute_thresholds(fit.ls$alert_L %>% filter(!grepl("Null", model)), 
                     0, 1, 0.05, byPrevAlert=T, cores=12),
  compute_thresholds(fit.ls$alert_L %>% filter(grepl("Null", model)), 
                     0, 1, 0.05, byPrevAlert=F, cores=12)
)
gc()
opt.J <- list(
  thresh.fit[[1]] %>% filter(!is.na(J)) %>%
    group_by(y, model, PCA, covSet, prevAlert) %>%
    arrange(desc(J)) %>% slice_head(n=1) %>% ungroup %>%
    select(y, model, PCA, covSet, thresh, J, prevAlert) %>%
    rename(optJ=thresh, J_fit=J),
  thresh.fit[[2]] %>% filter(!is.na(J)) %>%
    group_by(y, model, PCA, covSet) %>%
    arrange(desc(J)) %>% slice_head(n=1) %>% ungroup %>%
    select(y, model, PCA, covSet, thresh, J) %>%
    rename(optJ=thresh, J_fit=J)
)
opt.F1 <- list(
  thresh.fit[[1]] %>% filter(!is.na(F1)) %>%
    group_by(y, model, PCA, covSet, prevAlert) %>%
    arrange(desc(F1)) %>% slice_head(n=1) %>% ungroup %>%
    select(y, model, PCA, covSet, thresh, F1, prevAlert) %>%
    rename(optF1=thresh, F1_fit=F1),
  thresh.fit[[2]] %>% filter(!is.na(F1)) %>%
    group_by(y, model, PCA, covSet) %>%
    arrange(desc(F1)) %>% slice_head(n=1) %>% ungroup %>%
    select(y, model, PCA, covSet, thresh, F1) %>%
    rename(optF1=thresh, F1_fit=F1)
)
oos.ls$alert_L <- bind_rows(
  oos.ls$alert_L %>%
    filter(!grepl("Null", model)) %>%
    left_join(opt.J[[1]] %>% select(-J_fit)) %>%
    left_join(opt.F1[[1]] %>% select(-F1_fit)),
  oos.ls$alert_L %>%
    filter(grepl("Null", model)) %>%
    left_join(opt.J[[2]] %>% select(-J_fit)) %>%
    left_join(opt.F1[[2]] %>% select(-F1_fit))
  ) %>%
  mutate(predJ=factor(if_else(prA1 > optJ, "A1", "A0"), levels=c("A0", "A1")),
         predF=factor(if_else(prA1 > optF1, "A1", "A0"), levels=c("A0", "A1")))



# metrics -----------------------------------------------------------------

# AUC
oos.ls$alert_L %>% 
  filter(model != "perfect") %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  roc_auc(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(desc(.estimate)) %>% slice_head(n=1) %>%
  ungroup %>% arrange(desc(.estimate))
oos.ls$alert_L %>% 
  filter(model != "perfect") %>%
  filter(model != "Null (grand)") %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(desc(.estimate)) %>% slice_head(n=1) %>%
  ungroup %>% arrange(desc(.estimate))


# Scaled PR-AUC
oos.ls$alert_L %>% 
  filter(model != "perfect") %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(model) %>%
  mutate(propImprove=(.estimate-first(.estimate))/(1-first(.estimate))) %>%
  group_by(y) %>% arrange(desc(.estimate)) %>% slice_head(n=1) %>%
  ungroup %>% arrange(desc(.estimate))
oos.ls$alert_L %>% 
  filter(model %in% c("Null (grand)", "Ens-LogitWtMn")) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(model) %>%
  mutate(propImprove=(.estimate-first(.estimate))/(1-first(.estimate))) %>%
  ungroup %>% filter(model != "Null (grand)") %>%
  summary


oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  roc_auc(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0.5, 1) + labs(x="", y="ROC-AUC (oos)") + 
  facet_grid(.~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>% 
  filter(model != "perfect") %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(model) %>%
  mutate(.estimate=(.estimate-first(.estimate))/(1-first(.estimate))) %>%
  ungroup %>% 
  filter(model != "Null (grand)") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Scaled PR-AUC (oos)") + 
  facet_grid(.~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>% 
  # filter(grepl("^d", covSet)) %>%
  # filter(model %in% c("Ridge GLM", "MARS")) %>%
  # filter(model %in% c("Bayes-L1", "Bayes-NL1")) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="PR-AUC (oos)") + 
  facet_grid(.~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L %>% 
  na.omit %>%
  group_by(y, model, covSet, PCA, prevAlert) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="PR-AUC (oos)") + 
  facet_grid(prevAlert~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>% 
  filter(model != "perfect") %>%
  filter(model != "Null (grand)") %>%
  na.omit %>%
  group_by(y, model, covSet, PCA, prevAlert) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  group_by(y, prevAlert) %>% arrange(model) %>%
  mutate(.estimate=(.estimate-first(.estimate))/(1-first(.estimate))) %>%
  ungroup %>% 
  filter(model != "Null (grand)") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Scaled PR-AUC (oos)") + 
  facet_grid(prevAlert~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>% 
  group_by(y, model, PCA, covSet) %>%
  precision(predF, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Precision: P[y=1|y'=1] (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>% 
  group_by(y, model, PCA, covSet) %>%
  recall(predF, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Recall: P[y'=1|y=1] (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>% 
  group_by(y, model, PCA, covSet) %>%
  f_meas(predF, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="F1 (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>% 
  group_by(y, model, PCA, covSet) %>%
  j_index(predJ, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="J (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>%
  # filter(grepl("d", covSet)) %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  group_by(y, model, PCA, covSet) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
  group_by(y) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  filter(!is.na(R2)) %>%
  # filter(model %in% c("Ridge GLM", "MARS")) %>%
  # filter(model %in% c("Bayes-L1", "Bayes-NL1")) %>%
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
oos.ls$alert_L %>%
  # filter(grepl("d", covSet)) %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(PCA != "lnN") %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  group_by(y, model, PCA, covSet) %>%
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n())) %>%
  filter(!is.na(R2_mn)) %>%
  # filter(model %in% c("Ridge GLM", "MARS")) %>%
  # filter(model %in% c("Bayes-L1", "Bayes-NL1")) %>%
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
oos.ls$alert_L %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(PCA != "lnN") %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  ggplot(aes(covSet, R2, fill=model)) + 
  geom_boxplot() +
  scale_fill_manual(values=mod_cols) +
  labs(x="Model", y="McFaddens R2 mn, 2se among sites (oos)") + 
  ylim(-1, 1) + 
  facet_wrap(~y) +
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")




oos.ls$alert_L %>% 
  # filter(grepl("Neural", model)) %>%
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=PCA, group=PCA)) + geom_point() + geom_line() +
  facet_grid(y~model)
oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA, Spatial) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(PCA, .estimate, fill=Spatial)) + 
  geom_boxplot() +
  facet_wrap(~y, nrow=2)
oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA, EnvInt) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(PCA, .estimate, fill=EnvInt)) + 
  geom_boxplot() +
  facet_wrap(~y, nrow=2)
oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA, AutoInt) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(PCA, .estimate, fill=AutoInt)) + 
  geom_boxplot() +
  facet_wrap(~y, nrow=2)
oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA, EnvDelta) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(PCA, .estimate, fill=EnvDelta)) + 
  geom_boxplot() +
  facet_wrap(~y, nrow=2)


oos.ls$alert_L %>%
  filter(grepl("d", covSet)) %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  group_by(y, model, PCA, covSet) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
  group_by(y) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  filter(!is.na(R2)) %>%
  mutate(R2=if_else(R2 < 0, 0, R2)) %>%
  ggplot(aes(covSet, R2, colour=PCA, group=PCA)) + geom_point() + geom_line() +
  facet_grid(y~model)





oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, group=model)) + 
  geom_point() + geom_line() +
  scale_colour_manual(values=mod_cols) +
  facet_grid(PCA~y) + labs(x="", y="PR-AUC")

oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, PCA, covSet) %>%
  f_meas(predF, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, group=model)) + 
  geom_point() + geom_line() +
  scale_colour_manual(values=mod_cols) +
  facet_grid(PCA~y) + labs(x="", y="F1")

oos.ls$alert_L %>% 
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  group_by(y, model, PCA, covSet) %>%
  j_index(predJ, truth=alert, event_level="second") %>%
  ggplot(aes(covSet, .estimate, colour=model, group=model)) + 
  geom_point() + geom_line() +
  scale_colour_manual(values=mod_cols) +
  facet_grid(PCA~y) + labs(x="", y="J")

oos.ls$alert_L %>%
  filter(grepl("^d", covSet)) %>%
  na.omit %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  group_by(y, model, PCA, covSet) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
  group_by(y) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  filter(!is.na(R2)) %>%
  mutate(R2=if_else(R2 < 0, 0, R2)) %>%
  ggplot(aes(covSet, R2, colour=model, group=model)) + 
  geom_point() + geom_line() +
  scale_colour_manual(values=mod_cols) +
  facet_grid(PCA~y)


oos.ls$alert_L %>%
  filter(model != "perfect") %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  group_by(y, model, PCA, covSet) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
  group_by(y) %>%
  arrange(model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  arrange(desc(R2)) %>%
  slice_head(n=1) %>%
  ungroup %>%
  arrange(desc(R2))


oos.ls$alert_L %>%
  filter(model != "perfect") %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(PCA != "lnN") %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  group_by(y, model, PCA, covSet) %>%
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n()), nA1=sum(nA1), mnPA1=mean(pA1)) %>%
  ungroup %>%
  filter(!is.na(R2_mn)) %>%
  group_by(y) %>%
  arrange(R2_mn) %>% 
  slice_tail(n=1) %>%
  arrange(desc(R2_mn))


oos.ls$alert_L %>%
  filter(model != "perfect") %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(PCA != "lnN") %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  group_by(y, model, PCA, covSet) %>%
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n()), nA1=sum(nA1), mnPA1=mean(pA1)) %>%
  ungroup %>%
  filter(!is.na(R2_mn)) %>%
  ggplot(aes(mnPA1, R2_mn)) + 
  geom_point() + facet_wrap(~model)
siteR2 <- oos.ls$alert_L %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(PCA != "lnN") %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") 
siteR2 %>%
  filter(pA1 > 0 & pA1 < 1) %>%
  ggplot(aes(pA1, R2, colour=model)) + ylim(-1,1) + 
  scale_colour_manual(values=mod_cols) +
  geom_point(alpha=0.6, shape=1) + 
  stat_smooth(se=F, method="lm") + 
  facet_wrap(~y) +
  scale_x_continuous(trans="logit") +
  theme(legend.position="none")
siteR2 %>% 
  group_by(y, siteid) %>%
  summarise(mnR2=median(R2, na.rm=T), pA1=first(pA1)) %>%
  ggplot(aes(pA1, mnR2)) + 
  geom_point() +
  facet_wrap(~y)

siteAUC <- oos.ls$alert_L %>%
  filter(!grepl("grand", model)) %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(pA1=mean(alert=="A1")) %>%
  ungroup %>%
  filter(pA1 > 0) %>%
  group_by(y, model, PCA, covSet, siteid, pA1) %>%
  roc_auc(prA1, truth=alert, event_level="second") 
siteAUC %>%
  ggplot(aes(covSet, .estimate, fill=model)) + ylim(0,1) + 
  scale_fill_manual(values=mod_cols) +
  geom_boxplot() + 
  facet_wrap(~y) +
  theme(legend.position="none")
siteAUC %>%
  ggplot(aes(pA1, .estimate, colour=model)) + ylim(0,1) + 
  scale_colour_manual(values=mod_cols) +
  geom_point(alpha=0.6, shape=1) + 
  stat_smooth(se=F, method="lm") + 
  facet_wrap(~y) +
  scale_x_continuous(trans="logit") +
  theme(legend.position="none")
siteAUC %>%
  ggplot(aes(siteid, .estimate, colour=pA1, group=siteid)) +
  geom_boxplot() + 
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
siteAUC %>% 
  group_by(y, siteid) %>%
  summarise(mnAUC=median(.estimate, na.rm=T), pA1=first(pA1)) %>%
  ggplot(aes(pA1, mnAUC)) + 
  geom_point() +
  facet_wrap(~y)

sitePRAUC <- oos.ls$alert_L %>%
  filter(!grepl("grand", model)) %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(pA1=mean(alert=="A1")) %>%
  ungroup %>%
  filter(pA1 > 0) %>%
  group_by(y, model, PCA, covSet, siteid, pA1) %>%
  average_precision(prA1, truth=alert, event_level="second") 
sitePRAUC %>%
  ggplot(aes(covSet, .estimate, fill=model)) + ylim(0,1) + 
  scale_fill_manual(values=mod_cols) +
  geom_boxplot() + 
  facet_wrap(~y) +
  theme(legend.position="none")
sitePRAUC %>%
  ggplot(aes(pA1, .estimate, colour=model)) + ylim(0,1) + 
  scale_colour_manual(values=mod_cols) +
  geom_point(alpha=0.6, shape=1) + 
  stat_smooth(se=F, method="lm") + 
  facet_wrap(~y) +
  scale_x_continuous(trans="logit") +
  theme(legend.position="none")
sitePRAUC %>%
  ggplot(aes(siteid, .estimate, colour=pA1, group=siteid)) +
  geom_boxplot() + 
  scale_colour_viridis_c() +
  facet_wrap(~y, scales="free_x")
sitePRAUC %>% 
  group_by(y, siteid) %>%
  summarise(mnAUC=median(.estimate, na.rm=T), pA1=first(pA1)) %>%
  ggplot(aes(pA1, mnAUC)) + 
  geom_point() +
  facet_wrap(~y)



oos.ls$alert_L %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(PCA != "lnN") %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  mutate(anyA1=pA1 > 0) %>%
  ggplot(aes(anyA1, R2, fill=model)) + ylim(-1,1) + 
  scale_fill_manual(values=mod_cols) +
  geom_boxplot() +
  facet_wrap(~y) +
  theme(legend.position="none")



oos.ls$alert_L %>% 
  na.omit %>%
  # filter(model %in% c("Null (±2wk avg)", "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "perfect")) %>%
  group_by(y, model, PCA, covSet) %>% 
  pr_curve(prA1, truth=alert, event_level="second") %>% 
  filter(recall > 0) %>%
  ggplot(aes(recall, precision, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y) + ggtitle("oos")
oos.ls$alert_L %>% 
  filter(model %in% c("Null (±2wk avg)", "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "perfect")) %>%
  group_by(y, model, PCA, covSet) %>% 
  gain_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(.percent_tested, .percent_found, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)
oos.ls$alert_L %>% 
  # filter(model %in% c("Null (±2wk avg)", "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "perfect")) %>%
  na.omit %>%
  group_by(y, model, PCA, covSet) %>% 
  roc_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(1-specificity, sensitivity, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)

fit.ls$alert_L %>% 
  na.omit %>%
  # filter(model %in% c("Null (±2wk avg)", "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge")) %>%
  group_by(y, model, PCA, covSet) %>% 
  pr_curve(prA1, truth=alert, event_level="second") %>% 
  filter(recall > 0) %>%
  ggplot(aes(recall, precision, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y) + ggtitle("fit")
fit.ls$alert_L %>% 
  filter(model %in% c("Null (±2wk avg)", "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge")) %>%
  group_by(y, model, PCA, covSet) %>% 
  gain_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(.percent_tested, .percent_found, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)
fit.ls$alert_L %>% 
  filter(model %in% c("Null (±2wk avg)", "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge")) %>%
  na.omit %>%
  group_by(y, model, PCA, covSet) %>% 
  roc_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(1-specificity, sensitivity, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)




oos.ls$alert_L %>% 
  filter(model != "perfect") %>%
  filter(!grepl("grand", model)) %>%
  na.omit() %>%
  group_by(y, model, PCA, covSet) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  group_by(y) %>%
  mutate(rank=row_number()) %>%
  group_by(model, covSet) %>%
  summarise(mdRank=median(rank),
            mnRank=mean(rank)) %>%
  ungroup %>%
  arrange(mnRank) %>%
  print(n=36)



rank.df <-  bind_rows(oos.ls$alert_L %>% 
                        filter(model != "perfect") %>%
                        filter(!grepl("grand", model)) %>%
                        na.omit() %>%
                        group_by(y, model, PCA, covSet) %>%
                        average_precision(prA1, truth=alert, event_level="second") %>%
                        arrange(y, desc(.estimate)) %>%
                        group_by(y) %>%
                        mutate(rank=row_number()),
                      oos.ls$alert_L %>% 
                        filter(model != "perfect") %>%
                        filter(!grepl("grand", model)) %>%
                        na.omit() %>%
                        group_by(y, model, PCA, covSet) %>%
                        roc_auc(prA1, truth=alert, event_level="second") %>%
                        arrange(y, desc(.estimate)) %>%
                        group_by(y) %>%
                        mutate(rank=row_number())) %>%
  bind_rows(oos.ls$alert_L %>%
              filter(model != "perfect") %>%
              filter(!grepl("grand", model)) %>%
              na.omit() %>%
              group_by(y, model, PCA, covSet) %>%
              f_meas(predF, truth=alert, event_level="second") %>%
              arrange(y, desc(.estimate)) %>%
              group_by(y) %>%
              mutate(rank=row_number())) %>%
  bind_rows(oos.ls$alert_L %>%
              filter(model != "perfect") %>%
              filter(!grepl("grand", model)) %>%
              na.omit() %>%
              group_by(y, model, PCA, covSet) %>%
              j_index(predJ, truth=alert, event_level="second") %>%
              arrange(y, desc(.estimate)) %>%
              group_by(y) %>%
              mutate(rank=row_number())) %>%
  bind_rows(oos.ls$alert_L %>%
              filter(model != "perfect") %>%
              mutate(prA1=if_else(prA1==0, 1e-5, prA1),
                     prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
              group_by(y, model, PCA, covSet) %>%
              mutate(alert=as.numeric(alert=="A1")) %>%
              summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
              group_by(y) %>%
              arrange(model) %>%
              mutate(.estimate=1 - LL/first(LL)) %>%
              select(-LL) %>%
              filter(model != "Null (grand)") %>%
              na.omit() %>%
              arrange(desc(.estimate)) %>%
              mutate(rank=row_number(),
                     .metric="McFadden's R2",
                     .estimator="binary"))

rank.df %>% group_by(model) %>%
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) %>%
  arrange(meanRank)
rank.df %>% group_by(model, .metric) %>%
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) %>%
  arrange(.metric, meanRank) %>%
  filter(grepl("Ens", model)) %>%
  print(n=21)

rank.df %>% group_by(covSet) %>%
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) %>%
  arrange(meanRank)
rank.df %>% group_by(covSet, PCA) %>%
  summarise(medRank=median(rank), 
            meanRank=mean(rank),
            bestRank=min(rank),
            worstRank=max(rank)) %>%
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

rank.df %>% 
  filter(grepl("^d", covSet)) %>%
  ggplot(aes(.metric, rank, fill=covSet)) + 
  geom_boxplot() +
  scale_fill_viridis_d() + 
  labs(x="", y="Rank")


ggplot(rank.df, aes(.metric, rank, fill=model)) + 
  geom_boxplot() +
  scale_fill_manual(values=mod_cols) + 
  labs(x="", y="Rank") +
  facet_grid(y~.metric, scales="free_x")




# probabilities -----------------------------------------------------------

oos.ls$alert_L %>%
  filter(grepl("2wk|Ens", model)) %>%
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_boxplot(outlier.size=0.5, outlier.shape=1, outlier.alpha=0.5) +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_wrap(~y) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>%
  filter(grepl("2wk|Ens", model)) %>%
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_violin(scale="width") +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_wrap(~y) + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

oos.ls$alert_L %>%
  filter(grepl("2wk|Ens", model)) %>%
  filter(!grepl("grand", model)) %>%
  group_by(y, alert, model) %>%
  get_intervals(prA1, "qi") %>%
  ggplot(aes(med, xmin=L05, xmax=L95, model, colour=alert)) + 
  geom_linerange(aes(xmin=L25, xmax=L75), position=position_dodge(width=0.5), linewidth=2) + 
  geom_linerange(aes(xmin=L10, xmax=L90), position=position_dodge(width=0.5), linewidth=1) + 
  geom_pointrange(position=position_dodge(width=0.5)) + 
  scale_colour_manual("Truth", values=c("grey", "red3")) +
  facet_wrap(~y, nrow=2) +
  xlim(0,1)

oos.ls$alert_L %>%
  filter(grepl("2wk|Ens", model)) %>%
  filter(!grepl("grand", model)) %>%
  group_by(y, prevAlert, alert, model) %>%
  get_intervals(prA1, "qi") %>%
  ggplot(aes(med, xmin=L05, xmax=L95, model, colour=alert)) + 
  geom_linerange(aes(xmin=L25, xmax=L75), position=position_dodge(width=0.5), linewidth=2) + 
  geom_linerange(aes(xmin=L10, xmax=L90), position=position_dodge(width=0.5), linewidth=1) + 
  geom_linerange(position=position_dodge(width=0.5), linewidth=0.25) + 
  geom_point(position=position_dodge(width=0.5), size=3, shape=3) + 
  scale_colour_manual("Truth", values=c("grey", "red3")) +
  facet_grid(prevAlert~y) +
  xlim(0,1)

oos.ls$alert_L %>%
  filter(grepl("2wk|Ens", model)) %>%
  ggplot(aes(prA1, model, fill=alert)) +
  ggridges::geom_density_ridges(rel_min_height=0.01, scale=1, alpha=0.5) +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_grid(.~y) +
  xlim(0,1)
oos.ls$alert_L %>%
  filter(grepl("2wk|Ens", model)) %>%
  ggplot(aes(prA1, model, fill=alert)) +
  ggridges::geom_density_ridges(rel_min_height=0.01, scale=1, alpha=0.5) +
  scale_fill_manual("Truth", values=c("grey", "red3")) +
  facet_grid(prevAlert~y) +
  xlim(0,1)






# variable importance -----------------------------------------------------

enet.f <- dirf("out/model_fits/", "_Ridge.rds", recursive=T)
rf.f <- dirf("out/model_fits/", "_RF.rds", recursive=T)
nn.f <- dirf("out/model_fits/", "_NN.rds", recursive=T)

vi.df <- bind_rows(
  enet.f %>% 
    map_dfr(~readRDS(.x) %>% tidy %>% mutate(f=.x)) %>% 
    filter(term!="(Intercept)") %>%
    group_by(f) %>% mutate(Importance=abs(estimate)/max(abs(estimate))) %>% ungroup %>%
    rename(Variable=term) %>% select(-estimate) %>%
    mutate(model="Ridge"), 
  rf.f %>%
    map_dfr(~readRDS(.x) %>% extract_fit_engine() %>% vip::vi(scale=T) %>% mutate(f=.x)) %>%
    group_by(f) %>% mutate(Importance=Importance/max(Importance)) %>% ungroup %>%
    mutate(model="RF"),
  nn.f %>%
    map_dfr(~readRDS(.x) %>% extract_fit_engine() %>% vip::vi(scale=T) %>% mutate(f=.x)) %>%
    group_by(f) %>% mutate(Importance=abs(Importance)/max(abs(Importance))) %>% ungroup %>%
    mutate(model="NN")
) %>%
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels),
         y=str_split_fixed(str_split_fixed(f, "/", 4)[,4], "_", 3)[,1]) 

vi.df %>%
  ggplot(aes(Importance, Variable, fill=model)) + 
  geom_bar(stat="identity", colour="grey30", position="dodge") + 
  scale_fill_manual(values=mod_cols) +
  facet_grid(.~y)


