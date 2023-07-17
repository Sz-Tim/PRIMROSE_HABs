# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Manuscript figures and analysis



# setup -------------------------------------------------------------------

pkgs <- c("tidyverse", "lubridate", "glue", "tidymodels", "sf")
lapply(pkgs, library, character.only=T)
theme_set(theme_classic())
source("code/00_fn.R")

theme_ms <- theme_classic() + 
  theme(legend.title=element_text(size=7),
        legend.text=element_text(size=7),
        legend.background=element_blank(),
        strip.text=element_text(size=7),
        axis.title=element_text(size=9),
        axis.text=element_text(size=7),
        axis.line=element_line(colour="grey30", linewidth=0.25),
        strip.background=element_rect(fill=NA, colour="grey30", linewidth=0.5))

y_i <- bind_rows(read_csv("data/i_hab.csv") |> arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv") |> arrange(abbr) |> mutate(type="tox")) |>
  mutate(fig_short=factor(fig_short, 
                          levels=c("Alexandrium", "PSP", "Dinophysis", "DSP", 
                                   "K. mikimotoi", "Ps-nitz.", 
                                   "Ps-nitz. del.", "Ps-nitz. ser.", "ASP")),
         fig_long=factor(fig_long, 
                          levels=c("Alexandrium", "PSP", "Dinophysis", "DSP (OA/DTXs/PTXs)", 
                                   "Karenia mikimotoi", "Pseudo-nitzschia", 
                                   "Pseudo-nitzschia delicatissima", "Pseudo-nitzschia seriata", "ASP")))

site_hab.df <- readRDS("data/site_hab_df.rds")
site_tox.df <- readRDS("data/site_tox_df.rds")
coast <- st_read("data/northAtlantic_footprint.gpkg") |>
  st_crop(st_bbox(c(xmin=10e3, xmax=500e3, ymin=50e4, ymax=1280e3), crs=st_crs(27700)))

mod_i <- tibble(levels=c("nullGrand", "null4wk", "nullAuto", "perfect",
                         "ens", "ensLogitMn", "ensGLM", "ensGLM2",
                         "ensHB", "ensRF", "ensRF2", 
                         "Ridge", "MARS", "NN", 
                         "RF", "Boost",
                         "HBL1"),
                labels=c("Null (Int)", "Null (\u00B12wk)", "Null (auto)", "perfect", 
                         "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "Ensemble", 
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
              rep("grey40", 7),
              "#b2df8a", "#33a02c", "#ff7f00", 
              "#cab2d6", "#6a3d9a",
              "#1f78b4") |>
  setNames(mod_i$labels)

out.fit <- readRDS("out/clean/out_fit.rds")
out.oos <- readRDS("out/clean/out_oos.rds")
rank.oos <- readRDS("out/clean/rank_oos.rds")







# Fig 1 -------------------------------------------------------------------

site_all.df <- bind_rows(site_hab.df |> select(sin, lon, lat) |> mutate(type="HABs"),
                         site_tox.df |> select(sin, lon, lat) |> mutate(type="Biotoxins"))
fig1a <- ggplot(coast) + 
  geom_sf(fill="#a6cee3", colour="#a6cee3") + 
  geom_point(data=site_all.df, aes(lon, lat, colour=type), shape=1, size=0.7) + 
  scale_x_continuous("Longitude", breaks=c(-8, -4)) +
  scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
  scale_colour_manual("Monitoring", values=c("#33a02c", "#ff7f00")) +
  guides(colour=guide_legend(override.aes=list(size=1))) +
  theme_ms +
  theme(legend.position=c(0.21, 0.88),
        legend.key.width=grid::unit(0.2, "cm"),
        legend.key.height=grid::unit(0.3, "cm"))

obs.df <- bind_rows(readRDS("data/0_init/hab_obs.rds") |> 
                      select(y, date, obsid, siteid, sin, alert, lnN, tl),
                    readRDS("data/0_init/tox_obs.rds") |> 
                      select(y, date, obsid, siteid, sin, alert, lnN, tl)) |>
  filter(!y %in% c("AZP", "YTX", "Prli")) |>
  filter(date < "2023-01-01")
obs.wkly <- obs.df |>
  mutate(wk=week(date)-1,
         year=year(date),
         dateStd=ymd("2021-01-01") + wk*7) |>
  group_by(y, dateStd, year) |>
  summarise(prA=mean(alert=="A1"),
            N=n()) |>
  ungroup() |>
  mutate(lo=pmax(pmin(prA - 1.96*sqrt((prA*(1-prA))/N), 1), 0),
         hi=pmax(pmin(prA + 1.96*sqrt((prA*(1-prA))/N), 1), 0)) |>
  left_join(y_i, by=c("y"="abbr"))

fig1b <- obs.wkly |> 
  filter(N > 3) |>
  ggplot() + 
  geom_line(aes(dateStd, prA, group=year, colour=year), linewidth=0.25) +
  scale_y_continuous("Observed proportion of alerts per week", breaks=c(0, 0.25, 0.5, 0.75)) + 
  scale_x_date("Date", date_breaks="3 months", date_labels="%b") +
  scale_colour_viridis_c("Year", option="G", end=0.9, breaks=c(2016, 2019, 2022)) +
  # scale_colour_brewer(type="qual", palette=2) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), nrow=2) + 
  guides(colour=guide_colorbar(barwidth=grid::unit(0.25, "cm"), 
                               barheight=grid::unit(1.25, "cm"))) +
  theme_ms + 
  theme(panel.border=element_rect(colour="grey30", linewidth=0.2, fill=NA),
        legend.position=c(0.93, 0.25),
        legend.title=element_text(vjust=1),
        panel.spacing=unit(1, "mm"))

fig1 <- ggpubr::ggarrange(fig1a, fig1b, nrow=1, widths=c(1, 1.8), labels="AUTO")
ggsave("figs/pub/fig_1.png", fig1, width=190, height=95, units="mm", dpi=400)


obs.wkly |>
  filter(N > 3) |>
  group_by(y, dateStd) |>
  summarise(mnPrA=mean(prA),
            sdPrA=sd(prA),
            N=n()) |>
  group_by(y) |>
  arrange(desc(mnPrA)) |>
  slice_head(n=1)




site_sum.df <- obs.df |>
  left_join(y_i |> select(abbr, type, fig_long), by=c("y"="abbr")) |>
  group_by(y, fig_long, type, sin) |>
  summarise(prA=mean(alert=="A1"),
            N=n()) |>
  ungroup() |>
  mutate(type=if_else(type=="hab", "HABs", "Biotoxins")) |>
  left_join(site_all.df)
figS1 <- ggplot(coast) + 
  geom_sf(fill="grey70", colour="grey70") + 
  geom_point(data=site_sum.df, aes(lon, lat, colour=prA), shape=1) + 
  scale_x_continuous("Longitude", breaks=c(-8, -4)) +
  scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
  scale_colour_viridis_c("(warnings+alerts)/N\n", option="rocket") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm")) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
ggsave("figs/pub/fig_S1.png", figS1, width=190, height=140, units="mm", dpi=400)

site_month_sum.df <- obs.df |>
  left_join(y_i |> select(abbr, type, fig_long), by=c("y"="abbr")) |>
  mutate(month=month(date, label=T, abbr=F)) |>
  group_by(y, fig_long, type, sin, month) |>
  summarise(prA=mean(alert=="A1"),
            N=n()) |>
  ungroup() |>
  mutate(type=if_else(type=="hab", "HABs", "Biotoxins")) |>
  left_join(site_all.df)
for(i in 1:nlevels(site_month_sum.df$fig_long)) {
  ii <- levels(site_month_sum.df$fig_long)[i]
  figS1 <- ggplot(coast) + 
    geom_sf(fill="grey70", colour="grey70") + 
    geom_point(data=filter(site_month_sum.df, fig_long==ii), aes(lon, lat, colour=prA, group=sin), shape=1) + 
    scale_x_continuous("Longitude", breaks=c(-8, -4)) +
    scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
    scale_colour_viridis_c("(warnings+alerts)/N", option="rocket", limits=c(0,1)) +
    ggtitle(ii) +
    guides(colour=guide_colourbar(title.position="top", title.hjust=0.5)) +
    theme_ms +
    theme(legend.position="bottom", 
          legend.key.height=unit(1, "mm")) + 
    facet_wrap(~month, nrow=3)
  ggsave(glue("figs/pub/fig_S1{letters[i]}.png"), figS1, width=140, height=190, units="mm", dpi=400)
}

library(gganimate)
animS1 <- ggplot(coast) + 
  geom_sf(fill="grey70", colour="grey70") + 
  geom_point(data=site_month_sum.df, aes(lon, lat, colour=prA, group=sin), shape=1) + 
  transition_states(month, transition_length=0, state_length=1) +
  scale_x_continuous("Longitude", breaks=c(-8, -4)) +
  scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
  scale_colour_viridis_c("(warnings+alerts)/N\n", option="rocket") +
  ggtitle("{closest_state}") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm")) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
anim_save("figs/pub/anim_S1.gif", animS1, nframes=12, fps=2, 
          width=190, height=145, res=400, units="mm")




# Fig 2 -------------------------------------------------------------------

fig2 <- rank.oos  |>
  filter(! covSet %in% c("nullGrand", "nullAuto", "ens", "ensLogitMn", "ensGLM2", "ensHB")) |>
  filter(.metric %in% c("PR-AUC", "R2-VZ", "F1", "Schoener's D", "Accuracy (F1)")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("PR-AUC", "R2-VZ", "F1", 
                                 "Accuracy (F1)", "Schoener's D"),
                        labels=c("PR-AUC", "R['VZ']^2", "F[1]", 
                                 "Accuracy~(F[1])", "D[overlap]"))) |>
  mutate(model=factor(model, levels=c("Null (\u00B12wk avg)", "Ens-Ridge", "Ridge", 
                                      "MARS", "NN", "RF", "XGB", "HB"),
                      labels=c("Null (\u00B12wk)", "Ensemble", "Ridge", 
                               "MARS", "NN", "RF", "XGB", "HB"))) |>
  arrange(y, .metric, rank) |>
  group_by(y, .metric) |>
  mutate(rank=min_rank(rank)) |>
  ungroup() |>
  ggplot(aes(model, rank, fill=model)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.7, outlier.alpha=0.5, size=0.25) +
  scale_fill_manual("Model", values=mod_cols, guide="none") + 
  labs(x="", y="Rank") +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  theme_ms +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())

ggsave("figs/pub/Fig_2.png", fig2, width=140, height=90, units="mm", dpi=300)





# Fig 3 -------------------------------------------------------------------

fig3 <- rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM")) |>
  filter(.metric %in% c("PR-AUC", "R2-VZ", "F1", "Schoener's D", "Accuracy (F1)")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("PR-AUC", "R2-VZ", "F1", 
                                 "Accuracy (F1)", "Schoener's D"),
                        labels=c("PR-AUC", "R['VZ']^2", "F[1]", 
                                 "Accuracy~(F[1])", "D[overlap]"))) |>
  mutate(model=factor(model, levels=c("Null (int.)", "Null (\u00B12wk avg)", "Ens-Ridge"),
                      labels=c("Null (Int)", "Null (\u00B12wk)", "Ensemble"))) |>
  ggplot(aes(model, .estimate, colour=y, group=y)) + 
  geom_point(shape=1) + geom_line() +
  scale_colour_brewer(type="qual", palette=3) +
  scale_y_continuous("Value", breaks=c(0, 0.5, 1), limits=c(0,1)) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  guides(colour=guide_legend(nrow=2)) + 
  theme_ms +
  theme(legend.position=c(0.325, 0.925),
        legend.title=element_blank(),
        legend.key.height=unit(3, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid=element_blank())

ggsave("figs/pub/Fig_3.png", fig3, width=140, height=90, units="mm", dpi=300)






# Fig 4 -------------------------------------------------------------------

fig4 <- out.oos |>
  left_join(y_i, by=c("y"="abbr")) |>
  filter(grepl("2wk|Ens-HB", model)) |>
  mutate(model=factor(model, levels=c("Null (\u00B12wk avg)", "Ens-HB"),
                      labels=c("Null (\u00B12wk)", "Ensemble"))) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_boxplot(outlier.size=0.25, outlier.shape=1, outlier.alpha=0.25, size=0.25) +
  scale_fill_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_y_continuous("Forecasted Pr(y=A1)", breaks=c(0, 0.5, 1)) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), nrow=2) + 
  theme_ms +
  theme(legend.position=c(0.925, 0.225),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_4.png", fig4, width=140, height=120, units="mm", dpi=300)

fig4_alt <- out.oos |>
  left_join(y_i, by=c("y"="abbr")) |>
  filter(grepl("2wk|Ens-Ridge2", model)) |>
  mutate(model=factor(model, levels=c("Null (\u00B12wk avg)", "Ens-Ridge2"),
                      labels=c("Null (\u00B12wk)", "Ensemble"))) |>
  group_by(y, fig_long, obsid, alert) |>
  arrange(model) |>
  summarise(deltaPrA1=last(prA1)-first(prA1)) |>
  ungroup() |>
  mutate(fig_long=forcats::lvls_revalue(fig_long, str_wrap(levels(y_i$fig_long), 17))) |>
  mutate(fig_long=forcats::lvls_reorder(fig_long, 9:1)) |>
  ggplot(aes(deltaPrA1, fig_long, fill=alert)) +
  geom_vline(xintercept=0) + 
  geom_boxplot(outlier.size=0.25, outlier.shape=1, outlier.alpha=0.25, size=0.25) +
  scale_fill_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_x_continuous("Change in Pr(y=A1)\nEnsemble - Null (\u00B12wk avg)", 
                     breaks=c(-1, 0, 1), limits=c(-1, 1)) +
  scale_y_discrete(position="right") +
  theme_ms +
  theme(legend.position=c(0.15, 0.15),
        axis.title.y=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_4_alt.png", fig4_alt, width=90, height=100, units="mm", dpi=300)


out.oos |>
  left_join(y_i, by=c("y"="abbr")) |>
  filter(grepl("2wk|Ens-HB", model)) |>
  mutate(model=factor(model, levels=c("Null (\u00B12wk avg)", "Ens-HB"),
                      labels=c("Null (\u00B12wk)", "Ensemble")),
         prevAlert=factor(prevAlert, levels=c("A0", "A1"),
                          labels=paste("t-1:", 0:1))) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_boxplot(outlier.size=0.25, outlier.shape=1, outlier.alpha=0.25, size=0.25) +
  scale_fill_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_y_continuous("Forecasted Pr(y=A1)", breaks=c(0, 0.5, 1)) +
  facet_grid(prevAlert~fig_long, labeller=label_wrap_gen(17)) + 
  theme_ms +
  theme(legend.position=c(0.925, 0.875),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())

out.oos |>
  left_join(y_i, by=c("y"="abbr")) |>
  filter(grepl("2wk|Ens-Ridge2", model)) |>
  mutate(model=factor(model, levels=c("Null (\u00B12wk avg)", "Ens-Ridge2"),
                      labels=c("Null (\u00B12wk)", "Ensemble"))) |>
  group_by(y, fig_long, obsid, prevAlert, alert) |>
  arrange(model) |>
  summarise(deltaPrA1=last(prA1)-first(prA1)) |>
  ungroup() |>
  mutate(fig_long=forcats::lvls_revalue(fig_long, str_wrap(levels(y_i$fig_long), 17))) |>
  mutate(fig_long=forcats::lvls_reorder(fig_long, 9:1)) |>
  ggplot(aes(deltaPrA1, fig_long, fill=alert)) +
  geom_vline(xintercept=0) + 
  geom_boxplot(outlier.size=0.25, outlier.shape=1, outlier.alpha=0.25, size=0.25) +
  scale_fill_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_x_continuous("Change in Pr(y=A1)\nEnsemble - Null (\u00B12wk avg)", 
                     breaks=c(-1, 0, 1), limits=c(-1, 1)) +
  scale_y_discrete(position="right") +
  facet_grid(~prevAlert) +
  theme_ms +
  theme(legend.position=c(0.15, 0.15),
        axis.title.y=element_blank(),
        panel.grid=element_blank())

# Fig 5 -------------------------------------------------------------------

# Relative effects....

