# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Initial model analysis



# setup -------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels")
lapply(pkgs, library, character.only=T)
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))
source("code/00_fn.R")

y_i <- bind_rows(read_csv("data/i_hab.csv") %>% arrange(abbr) %>% mutate(type="hab"),
                 read_csv("data/i_tox.csv") %>% arrange(abbr) %>% mutate(type="tox"))

mod_i <- tibble(levels=c("nullGrand", "null4wk", 
                         "ens", "ensLasso", "perfect",
                         "Ridge GLM", "ENet", "MARS", "RF", "NN", 
                         "SVMl", "SVMr", "Boost",
                         paste0("HBL", 1:3), paste0("HBN", 1:3)),
                labels=c("Null (grand)", "Null (\u00B12wk avg)", 
                         "Ensemble-WtMn", "Ensemble-glm", "perfect",
                         "Ridge GLM", "ENet GLM", "MARS", "RandForest", "NeuralNetwork",
                         "SVM-linear", "SVM-radial", "XGBoost",
                         paste0("Bayes-L", 1:3), 
                         paste0("Bayes-NL", 1:3)))
mod_cols <- c(rep("grey", 2),
              "black", "black", "black",
              "#b2df8a", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
              "#fdbf6f", "#ff7f00", "#6a3d9a",
              rep("#a6cee3", 3), rep("#1f78b4", 3)) %>%
  setNames(mod_i$labels)

# d.ls <- dirf("data/0_init/", "data_baked.*test") %>% map(readRDS)
fit.ls <- dirf("out/compiled", "_fit.rds") %>%
  map(~readRDS(.x)) %>% list_transpose() %>% map(bind_rows)
oos.ls <- dirf("out/compiled", "_oos.rds") %>%
  map(~readRDS(.x)) %>% list_transpose() %>% map(bind_rows)
fit.ls$alert_L <- fit.ls$alert %>% 
  mutate(perfect_A1=if_else(alert=="A0", 1e-5, 1-1e-5)) %>%
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") %>% 
  mutate(model=str_split_fixed(run, "_", 3)[,1], 
         resp=str_split_fixed(run, "_", 3)[,2],
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=str_remove(model, "PCA.") %>% str_remove("d.\\.")) %>% 
  # filter(model != "nullGrand") %>%
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels))  
oos.ls$alert_L <- oos.ls$alert %>% 
  mutate(perfect_A1=if_else(alert=="A0", 1e-5, 1-1e-5)) %>%
  pivot_longer(ends_with("_A1"), names_to="run", values_to="prA1") %>% 
  mutate(model=str_split_fixed(run, "_", 3)[,1], 
         resp=str_split_fixed(run, "_", 3)[,2],
         PCA=grepl("PCA", model),
         covSet=str_split_fixed(model, "\\.", 2)[,1],
         model=str_remove(model, "PCA.") %>% str_remove("d.\\.")) %>% 
  # filter(model != "nullGrand") %>%
  mutate(model=factor(model, levels=mod_i$levels, labels=mod_i$labels))
  
meta <- map_dfr(dirf("out/model_fits", "[A-Za-z]_tune", recursive=T), 
                ~readRDS(.x) %>% collect_metrics() %>% 
                  mutate(f=.x,
                         model=str_sub(str_split_fixed(f, "_alert_", 2)[,2], 1, -10),
                         y=str_split_fixed(str_split_fixed(f, "meta/", 2)[,2], "_alert", 2)[,1],
                         covSet=str_split_fixed(f, "/", 5)[,3])) %>%
  mutate(PCA=grepl("PCA", model),
         model=str_remove(model, "_PCA"))
meta %>% 
  group_by(f, .metric) %>% arrange(mean) %>% slice_tail(n=1) %>%
  ggplot(aes(mean, model, colour=covSet, shape=PCA)) + 
  xlim(0, 1) + 
  scale_colour_brewer(type="qual", palette=2) + 
  scale_shape_manual(values=c(1, 19)) + 
  geom_point(alpha=0.6, size=3) + 
  facet_grid(y~.metric)

meta %>% 
  mutate(covSet=str_sub(covSet, 1, 1)) %>%
  group_by(f, .metric) %>% arrange(mean) %>% slice_tail(n=1) %>%
  ggplot(aes(covSet, mean, colour=PCA, group=paste(model, PCA), linetype=model, shape=model)) + 
  ylim(0, 1) + 
  geom_line() + geom_point() +
  facet_grid(.metric~y)

meta %>% 
  group_by(f, .metric) %>% arrange(mean) %>% slice_tail(n=1) %>% ungroup %>%
  select(y, model, covSet, PCA, .metric, mean) %>%
  pivot_wider(names_from=".metric", values_from="mean") %>%
  arrange(y, model, covSet, PCA)

meta %>% filter(grepl("ENet", f)) %>% mutate(mixture=factor(mixture)) %>% 
  filter(.metric=="roc_auc")%>%
  ggplot(aes(penalty, mean, colour=mixture)) + ylim(0.5, 1) +
  geom_line() + scale_colour_viridis_d(option="turbo") + facet_wrap(~f)
meta %>% filter(grepl("RF", f)) %>% mutate(min_n=factor(min_n)) %>% 
  filter(.metric=="roc_auc") %>%
  ggplot(aes(trees, mean, colour=min_n)) + ylim(0.5, 1) +
  geom_line() + scale_colour_viridis_d(option="turbo") + facet_wrap(~f)
meta %>% filter(grepl("NN", f)) %>% mutate(epochs=factor(epochs)) %>% 
  filter(.metric=="roc_auc") %>%
  ggplot(aes(penalty, mean, colour=epochs)) + ylim(0.5, 1) +
  geom_line() + scale_colour_viridis_d(option="turbo") + facet_grid(f~hidden_units)

meta %>% 
  filter(.metric=="pr_auc") %>%
  select(y, mean, model, covSet, PCA) %>%
  group_by(y) %>% arrange(desc(mean)) %>% slice_head(n=1)

meta %>% 
  filter(.metric=="roc_auc") %>%
  select(y, mean, model, covSet, PCA) %>%
  group_by(y, model) %>% arrange(desc(mean)) %>% slice_head(n=1)




# HSS pr threshold --------------------------------------------------------
thresh.HSS.fit <- map_dfr(seq(0, 1, by=0.01), 
                      ~fit.ls$alert_L %>%
                        mutate(thresh=.x,
                               predBloom=as.numeric(prA1 >= thresh)) %>%
                        group_by(y, model, resp, covSet, thresh) %>%
                        summarise(TP=sum(alert=="A1" & predBloom==1),
                                  TN=sum(alert=="A0" & predBloom==0),
                                  FP=sum(alert=="A0" & predBloom==1),
                                  FN=sum(alert=="A1" & predBloom==0),
                                  N=n()) %>%
                        mutate(expCorrRand=1/N * ((TP+FN)*(TP+FP) + (TN+FN)*(TN+FP)),
                               HSS=((TP+TN)-expCorrRand)/(N-expCorrRand)) %>%
                        filter(model!="Null (grand)"))
thresh.HSS.oos <- map_dfr(seq(0, 1, by=0.01), 
                          ~oos.ls$alert_L %>%
                            mutate(thresh=.x,
                                   predBloom=as.numeric(prA1 >= thresh)) %>%
                            group_by(y, model, resp, covSet, thresh) %>%
                            summarise(TP=sum(alert=="A1" & predBloom==1),
                                      TN=sum(alert=="A0" & predBloom==0),
                                      FP=sum(alert=="A0" & predBloom==1),
                                      FN=sum(alert=="A1" & predBloom==0),
                                      N=n()) %>%
                            mutate(expCorrRand=1/N * ((TP+FN)*(TP+FP) + (TN+FN)*(TN+FP)),
                                   HSS=((TP+TN)-expCorrRand)/(N-expCorrRand)) %>%
                            filter(model!="Null (grand)"))
opt.HSS <- thresh.HSS.fit %>%
  group_by(y, model, resp, covSet) %>%
  arrange(desc(HSS)) %>%
  slice_head(n=1) %>%
  ungroup %>%
  select(y, model, resp, covSet, thresh, HSS) %>%
  rename(opt_thresh=thresh, HSS_fit=HSS)
opt.HSS_oss <- thresh.HSS.oos %>%
  full_join(., opt.HSS) %>%
  filter(thresh==opt_thresh)
oos.ls$alert_L <- oos.ls$alert_L %>% 
  left_join(., opt.HSS)

opt.HSS_oss %>%
  ggplot(aes(HSS, model, colour=model, shape=covSet)) + 
  geom_point() + geom_rug(sides="b") +
  scale_colour_manual(values=mod_cols) +
  facet_grid(y~.) +
  scale_x_continuous(limits=c(-0.1, 1), breaks=seq(0,1,0.25)) +
  labs(x="HSS", y="", title="Out-of-sample HSS at fitted optimal p(bloom)") +
  theme(panel.grid.minor.y=element_blank(),
        legend.position="bottom")
thresh.HSS.oos %>%
  ggplot(aes(thresh, HSS, group=paste(model, resp, covSet), colour=model)) +
  geom_line() +
  scale_colour_manual(values=mod_cols) +
  facet_grid(covSet~y) +
  scale_y_continuous(limits=c(-0.1, 1), breaks=c(0, 0.5, 1)) +
  scale_x_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1)) +
  labs(x="p(bloom) threshold", y="HSS", title="Out-of-sample HSS") +
  theme(panel.grid.minor.y=element_blank(),
        legend.position="bottom")




# metrics -----------------------------------------------------------------

# AUC
oos.ls$alert_L %>% 
  na.omit %>%
  group_by(y, model, resp, covSet, PCA) %>%
  roc_auc(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(desc(.estimate)) %>% slice_head(n=1) %>%
  ungroup %>% arrange(desc(.estimate))
oos.ls$alert_L %>% 
  filter(model != "Null (grand)") %>%
  na.omit %>%
  group_by(y, model, resp, covSet, PCA) %>%
  pr_auc(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(desc(.estimate)) %>% slice_head(n=1) %>%
  ungroup %>% arrange(desc(.estimate))
oos.ls$alert_L %>% 
  filter(model != "Null (grand)") %>%
  na.omit %>%
  group_by(y, model, resp, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  group_by(y) %>% arrange(desc(.estimate)) %>% slice_head(n=1) %>%
  ungroup %>% arrange(desc(.estimate))

oos.ls$alert_L %>% 
  na.omit %>%
  group_by(y, model, resp, covSet, PCA) %>%
  average_precision(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Avg Precision (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L %>% 
  na.omit %>%
  group_by(y, model, resp, covSet, PCA) %>%
  pr_auc(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="PR-AUC (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L %>% 
  na.omit %>%
  group_by(y, model, resp, covSet, PCA) %>%
  roc_auc(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0.5, 1) + labs(x="", y="ROC-AUC (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L %>% 
  na.omit %>%
  group_by(y, model, resp, covSet, PCA) %>%
  gain_capture(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Gain-capture (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L %>%
  group_by(y, model, resp, covSet) %>%
  classification_cost(prA1, truth=alert, event_level="second", 
                                 costs=tibble(truth=c("A0", "A0", "A1", "A1"),
                                              estimate=c("A0", "A1", "A0", "A1"),
                                              cost=c(0, 1, 2, 0))) %>%
  arrange(y, desc(.estimate)) %>%
  ggplot(aes(model, .estimate, colour=model, shape=resp)) + 
  geom_point(size=2) + 
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  labs(x="", y="classification cost (oos)") + 
  facet_grid(.~y) +
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")
oos.ls$alert_L %>% 
  mutate(prA1=factor(prA1 > opt_thresh, levels=c("FALSE", "TRUE"), labels=c("A0", "A1"))) %>%
  # filter(resp != "lnN") %>%
  group_by(y, model, PCA, covSet) %>%
  kap(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(-1, 1) + labs(x="", y="kappa (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L %>% 
  mutate(prA1=factor(prA1 > opt_thresh, levels=c("FALSE", "TRUE"), labels=c("A0", "A1"))) %>%
  group_by(y, model, PCA, covSet) %>%
  f_meas(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  ggplot(aes(covSet, .estimate, colour=model, shape=PCA)) + 
  geom_point() + 
  geom_line(aes(group=paste(model, PCA))) +
  geom_rug(sides="l") +
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(19, 1)) +
  ylim(0, 1) + labs(x="", y="Fmeas (oos)") + 
  facet_wrap(~y) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
oos.ls$alert_L %>%
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
  arrange(y, desc(R2)) %>%
  ggplot(aes(covSet, R2, colour=model, linetype=PCA)) + geom_point() +
  geom_line(aes(group=paste(model, PCA)), linewidth=0.8) + 
  scale_colour_manual(values=mod_cols) +
  scale_shape_manual(values=c(2,1)) +
  ylim(0, NA) + labs(x="Model", y="McFaddens R2 (oos)") + 
  facet_grid(.~y) +
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")
oos.ls$alert_L %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(resp != "lnN") %>%
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
  ggplot(aes(covSet, R2_mn, colour=model)) + geom_point() + 
  geom_line(aes(linetype=PCA, group=paste(model, PCA)), linewidth=0.8) + 
  # geom_errorbar(aes(ymin=R2_mn-2*se, ymax=R2_mn+2*se), width=0.5) +
  scale_colour_manual(values=mod_cols) +
  labs(x="Model", y="McFaddens R2 mn, 2se among sites (oos)") + 
  facet_grid(.~y) + ylim(0, NA) + 
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")
oos.ls$alert_L %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(resp != "lnN") %>%
  group_by(y, model, PCA, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T))) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  ggplot(aes(model, R2, fill=model)) + 
  geom_boxplot() +
  scale_fill_manual(values=mod_cols) +
  labs(x="Model", y="McFaddens R2 mn, 2se among sites (oos)") + 
  ylim(-1, 1) + 
  facet_grid(.~y) +
  theme(panel.grid.minor=element_blank(), 
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        legend.position="bottom")


oos.ls$alert_L %>%
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
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(resp != "lnN") %>%
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
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(resp != "lnN") %>%
  group_by(y, model, resp, covSet, siteid) %>%
  mutate(alert=as.numeric(alert=="A1")) %>%
  summarise(LL=sum(dbinom(alert, 1, prA1, log=T)), 
            nA1=sum(alert==1),
            pA1=mean(alert)) %>%
  group_by(y, siteid) %>%
  arrange(y, model) %>%
  mutate(R2=1 - LL/first(LL)) %>%
  filter(model != "Null (grand)") %>%
  group_by(y, model, resp, covSet) %>%
  summarise(R2_mn=mean(R2), se=sd(R2)/sqrt(n()), nA1=sum(nA1), mnPA1=mean(pA1)) %>%
  ungroup %>%
  filter(!is.na(R2_mn)) %>%
  ggplot(aes(mnPA1, R2_mn)) + 
  geom_point() + facet_wrap(~model)
siteR2 <- oos.ls$alert_L %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(resp != "lnN") %>%
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


oos.ls$alert_L %>%
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1)) %>%
  # filter(resp != "lnN") %>%
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
  # filter(model!="Null (grand)") %>%
  group_by(y, model, resp, covSet) %>% 
  pr_curve(prA1, truth=alert, event_level="second") %>% 
  filter(recall > 0) %>%
  ggplot(aes(recall, precision, colour=model, group=paste(resp, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y) + ggtitle("oos")
oos.ls$alert_L %>% 
  group_by(y, model, resp, covSet) %>% 
  gain_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(.percent_tested, .percent_found, colour=model, group=paste(resp, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)
oos.ls$alert_L %>% 
  filter(resp != "lnN") %>%
  na.omit %>%
  group_by(y, model, resp, covSet) %>% 
  roc_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(1-specificity, sensitivity, colour=model, group=paste(resp, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)

fit.ls$alert_L %>% 
  na.omit %>%
  # filter(model!="Null (grand)") %>%
  group_by(y, model, resp, covSet) %>% 
  pr_curve(prA1, truth=alert, event_level="second") %>% 
  filter(recall > 0) %>%
  ggplot(aes(recall, precision, colour=model, group=paste(resp, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y) + ggtitle("fit")
fit.ls$alert_L %>% 
  group_by(y, model, resp, covSet) %>% 
  gain_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(.percent_tested, .percent_found, colour=model, group=paste(resp, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)
fit.ls$alert_L %>% 
  na.omit %>%
  group_by(y, model, resp, covSet) %>% 
  roc_curve(prA1, truth=alert, event_level="second") %>% 
  ggplot(aes(1-specificity, sensitivity, colour=model, group=paste(resp, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=mod_cols) + 
  facet_wrap(~y)




oos.ls$alert_L %>% 
  filter(!grepl("grand", model)) %>%
  na.omit() %>%
  group_by(y, model, PCA, covSet) %>%
  roc_auc(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  group_by(y) %>%
  mutate(rank=row_number()) %>%
  group_by(model, covSet) %>%
  summarise(mdRank=median(rank),
            mnRank=mean(rank)) %>%
  ungroup %>%
  arrange(mnRank) %>%
  print(n=28)


oos.ls$alert_L %>% 
  filter(!grepl("grand", model)) %>%
  na.omit() %>%
  group_by(y, model, PCA, covSet) %>%
  roc_auc(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  group_by(y) %>%
  mutate(rank=row_number()) %>% 
  ggplot(aes(model, rank, fill=model)) + 
  geom_boxplot() +
  scale_fill_manual(values=mod_cols, guide="none") + 
  labs(x="", y="ROC-AUC Rank within response")

oos.ls$alert_L %>% 
  filter(!grepl("grand", model)) %>%
  na.omit() %>%
  group_by(y, model, PCA, covSet) %>%
  roc_auc(prA1, truth=alert, event_level="second") %>%
  arrange(y, desc(.estimate)) %>%
  group_by(y) %>%
  mutate(rank=row_number()) %>% 
  ggplot(aes(model, rank)) + 
  geom_boxplot(aes(fill=model)) + 
  facet_wrap(~y) +
  scale_fill_manual(values=mod_cols, guide="none") + 
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))











# probabilities -----------------------------------------------------------

oos.ls$alert_L %>%
  filter(model != "Null (grand)") %>%
  ggplot(aes(alert, prA1, fill=model)) + 
  geom_boxplot(outlier.size=0.5, outlier.shape=1, outlier.alpha=0.5) +
  scale_fill_manual(values=mod_cols) +
  facet_wrap(~y)

oos.ls$alert_L %>%
  filter(model != "Null (grand)") %>%
  ggplot(aes(alert, prA1, fill=model)) + 
  geom_violin(scale="width") +
  scale_fill_manual(values=mod_cols) +
  facet_wrap(~y)









# variable importance -----------------------------------------------------

enet.f <- dirf("out/model_fits", "_ENet.rds", recursive=T)
rf.f <- dirf("out/model_fits", "_RF.rds", recursive=T)
nn.f <- dirf("out/model_fits", "_NN.rds", recursive=T)

vi.df <- bind_rows(
  enet.f %>% 
    map_dfr(~readRDS(.x) %>% tidy %>% mutate(f=.x)) %>% 
    filter(term!="(Intercept)") %>%
    group_by(f) %>% mutate(Importance=abs(estimate)/max(abs(estimate))) %>% ungroup %>%
    rename(Variable=term) %>% select(-estimate) %>%
    mutate(model="ENet"), 
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


