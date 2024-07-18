# PRIMROSE: Predicting Risk and Impact of Harmful Events on the Aquaculture Sector
# HAB Forecasting in the UK and Ireland
# Tim Szewczyk
# Manuscript figures and analysis



# setup -------------------------------------------------------------------

pkgs <- c("tidyverse", "glue", "tidymodels", "sf", "scico")
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
theme_talk <- theme_classic() + 
  theme(legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.background=element_blank(),
        strip.text=element_text(size=10),
        axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        axis.line=element_line(colour="grey30", linewidth=0.25),
        strip.background=element_rect(fill=NA, colour="grey30", linewidth=0.5))

y_i <- bind_rows(read_csv("data/i_hab.csv") |> arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv") |> arrange(abbr) |> mutate(type="tox")) |>
  mutate(fig_short=factor(fig_short, 
                          levels=c("A. sp.", "PST", "D. sp.", "DST", 
                                   "P. sp.", "P. del.", "P. ser.", "AST", 
                                   "K. mik.")),
         fig_long=factor(fig_long, 
                          levels=c("Alexandrium", "PST", "Dinophysis", "DST (OA/DTXs/PTXs)",
                                   "Karenia mikimotoi",
                                   "Pseudo-nitzschia", 
                                   "Pseudo-nitzschia delicatissima", "Pseudo-nitzschia seriata", "AST")),
         y=factor(abbr, levels=c("Alsp", "PSP", "Disp", "DSP", "Pssp", 
                              "Psde", "Psse", "ASP", "Kami")))

site_hab.df <- readRDS("data/site_hab_df.rds")
site_tox.df <- readRDS("data/site_tox_df.rds")
coast <- st_read("data/northAtlantic_footprint.gpkg") |>
  st_crop(st_bbox(c(xmin=50e3, xmax=500e3, ymin=50e4, ymax=1280e3), crs=st_crs(27700)))

UK_bbox <- st_multipoint(rbind(c(-8, 61.25),
                               c(-8, 54),
                               c(2, 54),
                               c(2, 61.25),
                               c(-8,61.25))) %>%
  st_cast("POLYGON") %>%
  st_sfc() %>%
  st_sf() %>%
  st_set_crs(4326) %>%
  st_transform(27700)

coast.path <- ifelse(Sys.info()["sysname"]=="Windows",
                     "../../00_gis/fromMTB/GIS/Layers/Physical/",
                     "../gis/coastlines/")
coast.sf <- st_read(paste0(coast.path, "simplified_land_polygons.shp")) |>
  st_transform(27700) |>
  st_make_valid() |>
  st_crop(UK_bbox) |>
  mutate(FID=1) |>
  st_union()

mod_i <- tibble(levels=c("nullGrand", "null4wk", "nullAuto", "perfect",
                         "ens", "ensLogitMn", "ensGLM", "ensGLM2",
                         "ensHB", "ensRF", "ensRF2", 
                         "HBL1", "Ridge", "MARS", "NN", 
                         "RF", "Boost"),
                labels=c("Null[0]", "Null[Date]", "Null[auto]", "perfect", 
                         "Ens-WtMn", "Ens-LogitWtMn", "Ens-Ridge", "Ensemble", 
                         "Ens-HB", "Ens-RF", "Ens-RF2", 
                         "HB", "Ridge", "MARS", "NN",
                         "RF", "XGB"))
mod_cols <- c(rep("grey", 3), "grey20",
              rep("grey40", 7),
              scico(6, palette="glasgow", begin=0.2, end=0.8, alpha=0.9)) |>
              # "#1f78b4", "#b2df8a", "#33a02c", "#ff7f00", 
              # "#cab2d6", "#6a3d9a") |>
  setNames(mod_i$labels)

out.oos <- readRDS("out/clean/out_oos.rds") |>
  mutate(y=factor(y, levels=levels(y_i$y)))
rank.oos <- readRDS("out/clean/rank_oos.rds") |>
  mutate(model=factor(model, levels=levels(model), labels=mod_i$labels),
         y=factor(y, levels=levels(y_i$y)))
rankPrev.oos <- readRDS("out/clean/rankPrevA_oos.rds") |>
  mutate(model=factor(model, levels=levels(model), labels=mod_i$labels),
         y=factor(y, levels=levels(y_i$y)))
vi_pt_df <- readRDS("out/clean/vi_df.rds") |>
  mutate(y=factor(y, levels=levels(y_i$y))) |>
  group_by(y, Variable, VariableType, varTypeClean) |>
  summarise(mnImp=mean(ImpWt, na.rm=T)) |>
  group_by(y) |>
  mutate(mnImp=mnImp/max(mnImp)*100,
         rank=min_rank(desc(mnImp)),
         top=(rank)/max(rank) <= 0.05) |>
  ungroup() |>
  left_join(y_i |> select(y, fig_short)) |>
  mutate(fig_short_rev=factor(fig_short, levels=rev(levels(fig_short))),
         varTypeClean=forcats::lvls_reorder(varTypeClean, c(1:4,6:9,5)))

obs.df <- bind_rows(readRDS("data/0_init/hab_obs.rds") |> 
                      select(y, date, obsid, siteid, sin, alert, lnN, tl),
                    readRDS("data/0_init/tox_obs.rds") |> 
                      select(y, date, obsid, siteid, sin, alert, lnN, tl)) |>
  filter(!y %in% c("AZP", "YTX", "Prli")) |>
  filter(date < "2023-01-01")






# Fig 1 -------------------------------------------------------------------

site_all.df <- bind_rows(read_csv("temp/fsa_sites_Scotland.csv") |> mutate(type="Phytoplankton"),
                         read_csv("temp/cefas_sites_Scotland.csv") |> mutate(type="Biotoxins")) |>
  filter(sin != "")
fig1a <- ggplot(coast) + 
  # geom_sf(fill="#a6cee3", colour="#a6cee3") + 
  geom_sf(fill=scico(15, palette="oslo")[12], colour=scico(15, palette="oslo")[11]) + 
  geom_point(data=site_all.df, aes(lon, lat, colour=type), shape=1, size=0.7) + 
  scale_x_continuous("Longitude", breaks=c(-6, -4, -2), expand=c(0, 0)) +
  scale_y_continuous("Latitude", breaks=c(55, 58, 61), expand=c(0, 0)) +
  scale_colour_manual("Monitoring", values=c("#33a02c", "#ff7f00")) +
  guides(colour=guide_legend(override.aes=list(size=1))) +
  theme_ms +
  theme(legend.position=c(0.23, 0.9),
        legend.key=element_rect(fill=NA, colour=NA),
        legend.key.width=grid::unit(0.2, "cm"),
        legend.key.height=grid::unit(0.3, "cm"))
Null_wkly <- out.oos |> 
  filter(model %in% c("Null[Date]", "Null[0]")) |> 
  mutate(model=factor(model, labels=c("italic(Null[0])", "italic(Null[Date])"))) |>
  mutate(yday=yday(date)) |> 
  group_by(yday, y, model) |> 
  summarise(prA=mean(prA1)) |> 
  ungroup() |>
  mutate(dateStd=ymd("2020-12-31") + yday,
         Null=model) |>
  rename(year=model) |>
  left_join(y_i) |> 
  select(dateStd, prA, year, fig_long, Null)
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
yr_cols <- c(`italic(Null[0])`="black",
             `italic(Null[Date])`="black",
             viridis::viridis(8) |> setNames(2015:2022))

fig1b <- obs.wkly |> 
  filter(N > 3) |>
  mutate(year=factor(year),
         Null="Obs") |>
  bind_rows(Null_wkly) |>
  ggplot() + 
  geom_line(aes(dateStd, prA, colour=year, linetype=Null, alpha=Null, linewidth=Null)) +
  scale_y_continuous("Weekly proportion of sites ≥ amber", breaks=c(0, 0.25, 0.5, 0.75)) + 
  scale_x_date(date_breaks="3 months", date_labels="%b") +
  # scale_x_date(breaks=seq(ymd("2021-01-01"), ymd("2021-12-01"), by="month"), 
  #              labels=str_sub(month.abb, 1, 1)) +
  scale_linetype_manual(values=c(3, 1, 1), guide="none") +
  scale_linewidth_manual(values=c(0.25, 0.4, 0.25), guide="none") +
  scale_alpha_manual(values=c(1, 1, 0.8), guide="none") +
  scale_colour_manual("Year", values=yr_cols, labels=parse_format()) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), nrow=2, scales="free_x") + 
  guides(colour=guide_legend(override.aes=list(linetype=c(rep(1, 8), 3, 1)))) +
  theme_ms + 
  theme(panel.border=element_rect(colour="grey30", linewidth=0.2, fill=NA),
        legend.position=c(0.915, 0.25),
        legend.key.height=unit(0.25, "cm"),
        legend.key.width=unit(0.45, "cm"),
        legend.title=element_blank(),
        panel.spacing=unit(1, "mm"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=6))

fig1 <- ggpubr::ggarrange(fig1a, fig1b, nrow=1, widths=c(1, 2), labels="AUTO")
ggsave("figs/pub/Fig_1.png", fig1, width=190, height=95, units="mm", dpi=400)


obs.wkly |>
  filter(N > 3) |>
  group_by(y, dateStd) |>
  summarise(mnPrA=mean(prA),
            sdPrA=sd(prA),
            q10PrA=quantile(prA, probs=0.1),
            q90PrA=quantile(prA, probs=0.9),
            iqrPrA=IQR(prA),
            rangePrA=diff(range(prA)),
            N=n()) |>
  group_by(y) |>
  arrange(desc(mnPrA)) |>
  slice_head(n=1) |>
  ungroup() |>
  summary()


obs.mo.Site <- obs.df |>
  mutate(month=month(date),
         year=year(date),
         dateStd=ymd(paste0("2021-", month, "-01"))) |>
  group_by(y, dateStd, siteid) |>
  summarise(prA=mean(alert=="A1"),
            N=n()) |>
  ungroup() |>
  mutate(lo=pmax(pmin(prA - 1.96*sqrt((prA*(1-prA))/N), 1), 0),
         hi=pmax(pmin(prA + 1.96*sqrt((prA*(1-prA))/N), 1), 0)) |>
  left_join(y_i, by=c("y"="abbr"))
ggplot(obs.mo.Site, aes(dateStd, prA, group=siteid)) + geom_line(alpha=0.3) + 
  facet_wrap(~fig_long)
ggplot(obs.mo.Site, aes(dateStd, prA, group=dateStd)) + geom_boxplot() + 
  facet_wrap(~fig_long)


# Fig 3 -------------------------------------------------------------------

fig3_mod_labs <- c("italic(Null[Date])", "italic(Ensemble)")
fig3_talk_labs <- c("italic(Null[0])", "italic(Null[Date])", "italic(Ens.)")

fig3_metrics_df <- rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  mutate(model=factor(model, levels=c("Null[0]", "Null[Date]", "Ensemble"))) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="D[overlap]", -first(.estimate), (1-first(.estimate)))) |>
  ungroup()
fig3_df <- bind_rows(fig3_metrics_df, 
                     fig3_metrics_df |>
                       group_by(y, model) |>
                       summarise(across(where(is.numeric), mean)) |>
                       mutate(.metric="Mean") |>
                       ungroup()
) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]", "Mean"))) |>
  filter(!grepl("0", model)) |> 
  left_join(y_i |> select(y, fig_short, fig_long))
fig3_means <- fig3_df |>
  group_by(model, .metric) |>
  summarise(across(where(is.numeric), mean))
fig3 <- fig3_df |>
  ggplot(aes(model, std_gain, colour=fig_short)) + 
  geom_hline(yintercept=0, linewidth=0.25, colour="grey30") +
  geom_point() + geom_line(aes(group=fig_short)) + 
  geom_point(data=fig3_means, shape=1, size=3, colour="black") +
  scale_colour_brewer(type="qual", palette=3) +
  scale_x_discrete(labels=parse(text=fig3_talk_labs[-1])) +
  scale_y_continuous(expression(paste("Skill score (0: Null" [0]~~"1: Perfect)")),
                     breaks=c(0, 0.5, 1), limits=c(-0.05,1)) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  guides(colour=guide_legend(nrow=2)) + 
  theme_ms +
  theme(legend.position=c(0.45, 0.925),
        legend.title=element_blank(),
        legend.background=element_rect(fill="white", colour=NA),
        legend.key.height=unit(3, "mm"),
        legend.key.width=unit(4, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.2),
        panel.grid.minor.y=element_line(colour="grey90", linewidth=0.2))
ggsave("figs/pub/Fig_3.png", fig3, width=140, height=90, units="mm", dpi=300)








fig3 <- rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  mutate(model=factor(model, levels=c("Null[0]", "Null[Date]", "Ensemble"))) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="D[overlap]", -first(.estimate), (1-first(.estimate)))) |>
  ungroup() |>
  filter(!grepl("Grand", covSet)) |>
  left_join(y_i |> select(y, fig_short)) |>
  ggplot(aes(model, .estimate, colour=fig_short, group=fig_short)) + 
  # geom_hline(yintercept=0, linewidth=0.25, colour="grey30") +
  geom_point(shape=1) + geom_line() + 
  scale_colour_brewer(type="qual", palette=3) +
  scale_x_discrete(labels=parse(text=fig3_mod_labs)) +
  ylab("Score") +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed", scales="free_y") + 
  guides(colour=guide_legend(nrow=2)) + 
  theme_ms +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        legend.key.height=unit(3, "mm"),
        legend.key.width=unit(4, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_3_scores.png", fig3, width=140, height=90, units="mm", dpi=300)

fig3 <- fig3_df |>
  ggplot(aes(model, std_gain, colour=fig_short)) + 
  geom_hline(yintercept=0, linewidth=0.25, colour="grey30") +
  geom_point() + geom_line(aes(group=fig_short)) + 
  geom_point(data=fig3_means, shape=1, size=4, colour="black") +
  scale_colour_brewer(type="qual", palette=3) +
  scale_x_discrete(labels=parse(text=fig3_talk_labs[-1])) +
  scale_y_continuous(expression(paste("Skill score (0: Null" [0]~~"1: Perfect)")),
                     breaks=c(0, 0.5, 1), limits=c(-0.05,1)) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  theme_talk +
  theme(legend.position="none", 
        legend.title=element_blank(),
        legend.key.height=unit(3, "mm"),
        legend.key.width=unit(4, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey80", linewidth=0.25),
        panel.grid.minor.y=element_line(colour="grey80", linewidth=0.25))
ggsave("figs/pub/Fig_3_talk.png", fig3, width=160, height=120, units="mm", dpi=300)

fig3 <- fig3_df |>
  ggplot(aes(model, std_gain, colour=fig_short)) + 
  geom_hline(yintercept=0, linewidth=0.25, colour="grey30") +
  geom_point() + geom_line(aes(group=fig_short)) + 
  geom_point(data=fig3_means, shape=1, size=3, colour="black") +
  scale_colour_brewer(type="qual", palette=3) +
  scale_x_discrete(labels=parse(text=fig3_talk_labs[-1])) +
  scale_y_continuous(expression(paste("Skill score (0: Null" [0]~~"1: Perfect)")),
                     breaks=c(0, 0.5, 1), limits=c(-0.05,1)) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  theme_talk +
  theme(legend.title=element_blank(),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey80"),
        panel.grid.minor.y=element_line(colour="grey80"))
ggsave("figs/pub/Fig_3_talk_legend.png", fig3, width=140, height=120, units="mm", dpi=300)

fig3_alt <- rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  mutate(model=factor(model, levels=c("Null[0]", "Null[Date]", "Ensemble"))) |>
  left_join(y_i |> select(y, fig_short)) |>
  ggplot(aes(model, .estimate, colour=fig_short, group=fig_short)) + 
  geom_point(shape=1) + geom_line() +
  scale_colour_brewer(type="qual", palette=3) +
  scale_y_continuous("Value", breaks=c(0, 0.5, 1), limits=c(0,1)) +
  scale_x_discrete(labels=parse(text=fig3_mod_labs)) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  guides(colour=guide_legend(nrow=2)) + 
  theme_ms +
  theme(legend.position=c(0.315, 0.925),
        legend.title=element_blank(),
        legend.key.height=unit(3, "mm"),
        legend.key.width=unit(4, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid=element_blank())

ggsave("figs/pub/Fig_3_alt.png", fig3_alt, width=140, height=90, units="mm", dpi=300)


sp_rankOrder <- rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(y, model) |>
  summarise(across(c("diff_v_null", "std_gain"), mean)) |>
  ungroup() |>
  filter(grepl("Ens", model)) |>
  arrange(desc(std_gain))

fig3_alt <- rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  mutate(model=factor(model, levels=c("Null[0]", "Null[Date]", "Ensemble"))) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="D[overlap]", -first(.estimate), (1-first(.estimate)))) |>
  filter(!grepl("0", model)) |>
  mutate(y=factor(y, levels=sp_rankOrder$y)) |>
  ggplot(aes(y, std_gain, colour=model, group=model)) + 
  geom_hline(yintercept=0, linewidth=0.25, colour="grey30") +
  geom_point(shape=1) + geom_line() + 
  scale_colour_manual(values=mod_cols) +
  scale_y_continuous("Skill score\n(0: Null; 1: Perfect)", 
                     breaks=c(0, 0.5, 1), limits=c(-0.05,1)) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  guides(colour=guide_legend(nrow=2)) + 
  theme_ms +
  theme(legend.position=c(0.325, 0.925),
        legend.title=element_blank(),
        legend.key.height=unit(3, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid=element_blank())

ggsave("figs/pub/Fig_3_alt2.png", fig3_alt, width=140, height=90, units="mm", dpi=300)

fig3_alt <- rankPrev.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  mutate(model=factor(model, levels=c("Null[0]", "Null[Date]", "Ensemble"))) |>
  group_by(.metric, y, prevAlert) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="D[overlap]", -first(.estimate), (1-first(.estimate)))) |>
  ggplot(aes(model, std_gain, colour=y, group=y)) + 
  geom_point(shape=1) + geom_line() +
  scale_colour_brewer(type="qual", palette=3) +
  scale_y_continuous("Value", breaks=c(0, 0.5, 1), limits=c(-0.1,1)) +
  scale_x_discrete(labels=parse(text=fig3_mod_labs)) +
  facet_grid(prevAlert~.metric, labeller="label_parsed") + 
  guides(colour=guide_legend(nrow=2)) + 
  theme_ms +
  theme(legend.position=c(0.325, 0.925),
        legend.title=element_blank(),
        legend.key.height=unit(3, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid=element_blank())

ggsave("figs/pub/Fig_3_alt3.png", fig3_alt, width=140, height=90, units="mm", dpi=300)










# Fig 4 -------------------------------------------------------------------

fig4_mod_labs <- c("italic(Null[Date])", "italic(Ensemble)")
fig4_talk_labs <- c("italic(Null[Date])", "italic(Ens.)")

fig4 <- out.oos |>
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  # filter(grepl("Date|Ens", model)) |>
  mutate(model=factor(model, levels=c("Null[Date]", "Ensemble"),
                      labels=c("Null[Date]", "Ensemble")),
         alert=factor(alert, levels=c("A0", "A1"), labels=c("< amber", "≥ amber"))) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_boxplot(outlier.size=0.25, outlier.shape=1, outlier.alpha=0.25, size=0.25) +
  scale_fill_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_y_continuous("Forecasted Pr(≥ amber)", breaks=c(0, 0.5, 1)) +
  scale_x_discrete(labels=parse(text=fig4_mod_labs)) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), nrow=2) + 
  theme_ms +
  theme(legend.position=c(0.925, 0.225),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_4.png", fig4, width=140, height=120, units="mm", dpi=300)

fig4 <- out.oos |>
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  # filter(grepl("Date|Ens", model)) |>
  mutate(model=factor(model, levels=c("Null[Date]", "Ensemble"),
                      labels=c("Null[Date]", "Ensemble")),
         alert=factor(alert, levels=c("A0", "A1"), labels=c("No alert", "Alert"))) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_boxplot(outlier.size=0.2, outlier.shape=1, outlier.alpha=0.1, size=0.25) +
  scale_fill_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_y_continuous("Forecasted Pr(Alert)", breaks=c(0, 0.5, 1)) +
  scale_x_discrete(labels=parse(text=fig4_talk_labs)) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), nrow=2) + 
  theme_talk +
  theme(legend.position=c(0.9, 0.125),
        # legend.title=element_blank(),
        legend.key.height=unit(5, "mm"),
        legend.key.width=unit(5, "mm"),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5), 
        axis.title.x=element_blank(),
        panel.grid=element_blank(), 
        panel.border=element_rect(colour="grey30", fill=NA))
ggsave("figs/pub/Fig_4_talk.png", fig4, width=170, height=120, units="mm", dpi=300)


fig4_alt <- out.oos |>
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  mutate(model=factor(model, levels=c("Null[Date]", "Ensemble"),
                      labels=c("Null[Date]", "Ensemble"))) |>
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
  scale_x_continuous("Change in Pr(y=A1)\nEnsemble - Null[date]", 
                     breaks=c(-1, 0, 1), limits=c(-1, 1)) +
  scale_y_discrete(position="right") +
  theme_ms +
  theme(legend.position=c(0.15, 0.15),
        axis.title.y=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_4_alt.png", fig4_alt, width=90, height=100, units="mm", dpi=300)


out.oos |>
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  mutate(prevAlert=factor(prevAlert, levels=c("A0", "A1"),
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
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
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
  scale_x_continuous("Change in Pr(y=A1)\nEnsemble - Null[date]", 
                     breaks=c(-1, 0, 1), limits=c(-1, 1)) +
  scale_y_discrete(position="right") +
  facet_grid(~prevAlert) +
  theme_ms +
  theme(legend.position=c(0.15, 0.15),
        axis.title.y=element_blank(),
        panel.grid=element_blank())

fig4_alt <- out.oos |>
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  group_by(fig_long, alert, model) |>
  get_intervals(prA1, "qi") |>
  ungroup() |>
  mutate(model=factor(model, levels=c("Null[Date]", "Ensemble"),
                      labels=c("Null[Date]", "Ensemble")),
         alert=factor(alert, levels=c("A0", "A1"), labels=c("No alert", "Alert"))) |>
  ggplot(aes(model, ymin=L025, ymax=L975, colour=alert)) + 
  geom_linerange(aes(ymin=L05, ymax=L95), position=position_dodge(width=0.5), linewidth=0.75) + 
  geom_linerange(aes(ymin=L10, ymax=L90), position=position_dodge(width=0.5), linewidth=1) + 
  geom_linerange(aes(ymin=L15, ymax=L85), position=position_dodge(width=0.5), linewidth=1.25) + 
  geom_linerange(aes(ymin=L20, ymax=L80), position=position_dodge(width=0.5), linewidth=1.5) + 
  geom_linerange(aes(ymin=L25, ymax=L75), position=position_dodge(width=0.5), linewidth=1.75) + 
  geom_linerange(aes(ymin=L30, ymax=L70), position=position_dodge(width=0.5), linewidth=2) + 
  geom_linerange(aes(ymin=L35, ymax=L65), position=position_dodge(width=0.5), linewidth=2.25) + 
  geom_linerange(aes(ymin=L40, ymax=L60), position=position_dodge(width=0.5), linewidth=2.5) + 
  geom_linerange(aes(ymin=L45, ymax=L55), position=position_dodge(width=0.5), linewidth=2.75) + 
  geom_linerange(position=position_dodge(width=0.5), linewidth=0.5) + 
  scale_colour_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_y_continuous("Forecasted Pr(Alert)", breaks=c(0, 0.5, 1)) +
  scale_x_discrete(labels=parse(text=fig4_mod_labs)) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), ncol=2) +
  theme_ms + 
  theme(legend.position=c(0.8, 0.05),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_4_alt2.png", fig4_alt, width=90, height=210, units="mm", dpi=300)

fig4_alt <- out.oos |>
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  group_by(fig_long, alert, model) |>
  get_intervals(prA1, "hdci") |>
  ungroup() |>
  mutate(model=factor(model, levels=c("Null[Date]", "Ensemble"),
                      labels=c("Null[Date]", "Ensemble")),
         alert=factor(alert, levels=c("A0", "A1"), labels=c("No alert", "Alert"))) |>
  ggplot(aes(model, ymin=L025, ymax=L975, colour=alert)) + 
  geom_linerange(aes(ymin=L05, ymax=L95), position=position_dodge(width=0.5), linewidth=0.75) + 
  geom_linerange(aes(ymin=L10, ymax=L90), position=position_dodge(width=0.5), linewidth=1) + 
  geom_linerange(aes(ymin=L15, ymax=L85), position=position_dodge(width=0.5), linewidth=1.25) + 
  geom_linerange(aes(ymin=L20, ymax=L80), position=position_dodge(width=0.5), linewidth=1.5) + 
  geom_linerange(aes(ymin=L25, ymax=L75), position=position_dodge(width=0.5), linewidth=1.75) + 
  geom_linerange(aes(ymin=L30, ymax=L70), position=position_dodge(width=0.5), linewidth=2) + 
  geom_linerange(aes(ymin=L35, ymax=L65), position=position_dodge(width=0.5), linewidth=2.25) + 
  geom_linerange(aes(ymin=L40, ymax=L60), position=position_dodge(width=0.5), linewidth=2.5) + 
  geom_linerange(aes(ymin=L45, ymax=L55), position=position_dodge(width=0.5), linewidth=2.75) + 
  geom_linerange(position=position_dodge(width=0.5), linewidth=0.5) + 
  scale_colour_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_y_continuous("Forecasted Pr(Alert)", breaks=c(0, 0.5, 1)) +
  scale_x_discrete(labels=parse(text=fig4_mod_labs)) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), ncol=2) +
  theme_ms + 
  theme(legend.position=c(0.8, 0.05),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_4_alt3.png", fig4_alt, width=90, height=210, units="mm", dpi=300)

fig4_alt <- out.oos |>
  left_join(y_i |> select(y, fig_long)) |>
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  # filter(grepl("Date|Ens", model)) |>
  mutate(model=factor(model, levels=c("Null[Date]", "Ensemble"),
                      labels=c("Null[Date]", "Ensemble")),
         alert=factor(alert, levels=c("A0", "A1"), labels=c("No alert", "Alert"))) |>
  ggplot(aes(model, prA1, fill=alert)) + 
  geom_violin(linewidth=0.25, scale="width", draw_quantiles=c(0.1, 0.5, 0.9)) +
  scale_fill_manual("Truth\n(out-of-sample)", values=c("grey", "red3")) +
  scale_y_continuous("Forecasted Pr(Alert)", breaks=c(0, 0.5, 1)) +
  scale_x_discrete(labels=parse(text=fig4_mod_labs)) +
  facet_wrap(~fig_long, labeller=label_wrap_gen(17), ncol=2) + 
  theme_ms +
  theme(legend.position=c(0.8, 0.05),
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())
ggsave("figs/pub/Fig_4_alt4.png", fig4_alt, width=90, height=210, units="mm", dpi=300)







# Fig 5 -------------------------------------------------------------------

library(GGally)

ggally_points_lims <- function(data, mapping, 
                               xlims=c(0, 1), ylims=c(0, 1)) {
  if(as.character(mapping$x)[2] == "MCC") {
    xlims <- c(-1, 1)
  }
  if(as.character(mapping$y)[2] == "MCC") {
    ylims <- c(-1, 1)
  }
  p <- ggplot(data=data, mapping=mapping) +
    geom_point(data=data |> filter(model=="Candidate"), 
               shape=1, size=0.3, alpha=0.4) +
    geom_point(data=data |> filter(model=="Null[Date]"), aes(fill=fig_short),
               colour="grey30", shape=21, size=1, alpha=0.8) +
    geom_point(data=data |> filter(model=="Ensemble"), aes(fill=fig_short),
               colour="grey30", shape=23, size=1.5, alpha=0.8) +
    {if(grepl("ROC", as.character(mapping$x)[2]) & 
        grepl("PR", as.character(mapping$y)[2])) { 
      geom_text(data=data |> group_by(model) |> summarise(x=0.15) |> ungroup() |>
                  mutate(y=1-row_number()*0.1, fig_short="DST"), 
                aes(label=model, x=x, y=y), parse=T,
                size=2, alpha=1, hjust=0, colour="grey20")
    }} +
    {if(grepl("ROC", as.character(mapping$x)[2]) &
        grepl("PR", as.character(mapping$y)[2])) {
      geom_point(data=data |> group_by(model) |> summarise(x=0.1) |> ungroup() |>
                   mutate(y=1-row_number()*0.1, fig_short="DST"),
                 aes(x=x, y=y, shape=model, size=model), alpha=1, 
                 colour="grey20", fill="grey50")
    }} +
    scale_shape_manual(values=c(21, 23, 1)) +
    scale_size_manual(values=c(1, 1, 0.5)) +
    scale_x_continuous(limits=xlims, 
                       breaks=seq(xlims[1], xlims[2], length.out=3), 
                       labels=seq(xlims[1], xlims[2], length.out=3)) +
    scale_y_continuous(limits=ylims, 
                       breaks=seq(ylims[1], ylims[2], length.out=3),
                       labels=seq(ylims[1], ylims[2], length.out=3)) +
    theme(panel.grid.major=element_line(colour="grey90", linewidth=0.2))
  return(p)
}


densityDiag_lims <- function(data, mapping, alpha=0.1, adjust=2, size=0.25, 
                             xlims=c(0, 1)) {
  if(as.character(mapping$x)[2] == "MCC") {
    xlims <- c(-1, 1)
  }
  p <- ggally_densityDiag(data, mapping, alpha=alpha, adjust=adjust, size=size) +
    scale_x_continuous(limits=xlims, 
                       breaks=seq(xlims[1], xlims[2], length.out=3), 
                       labels=seq(xlims[1], xlims[2], length.out=3)) +
    {if(grepl("ROC", as.character(mapping$x)[2])) { 
      geom_text(data=data |> group_by(fig_short) |> slice_head(n=1) |> ungroup() |>
                  mutate(y=13 - row_number()*1.25), 
                aes(label=fig_short, y=y), 
                x=0.15, size=2, alpha=1, hjust=0, colour="grey20")
    }} +
    {if(grepl("ROC", as.character(mapping$x)[2])) {
      geom_point(data=data |> group_by(fig_short) |> slice_head(n=1) |> ungroup() |>
                   mutate(y=13 - row_number()*1.25),
                 aes(y=y), x=0.1, size=1, alpha=1)
    }} +
    theme(panel.grid.major.x=element_line(colour="grey90", linewidth=0.2),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank())
  return(p)
}

theme_set(theme_classic() + 
            theme(strip.text=element_text(size=7),
                  axis.title=element_text(size=9),
                  axis.text=element_text(size=7),
                  axis.line=element_line(colour="grey30", linewidth=0.25),
                  strip.background=element_rect(fill=NA, colour="grey30", linewidth=0.5)))
fig5 <- rank.oos |> 
  filter(.metric %in% c("MCC", "R2-VZ", "PR-AUC", "ROC-AUC", "Schoener's D"),
         covSet %in% c("ensGLM2", "null4wk", paste0("d", 1:16))) |> 
  ungroup() |> 
  left_join(y_i |> select(y, fig_short)) |>
  select(fig_short, model, covSet, PCA, .metric, .estimate) |> 
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  arrange(.metric, fig_short, model) |>
  pivot_wider(names_from=".metric", values_from=".estimate") |> 
  select(-covSet, -PCA) |> 
  mutate(model=fct_collapse(model, 
                            Ensemble="Ensemble", 
                            `Null[Date]`="Null[Date]", 
                            other_level="Candidate")) |>
  ggpairs(columns=3:7, aes(colour=fig_short, fill=fig_short, group=fig_short,
                           size=model, alpha=model), 
          diag=list(
            continuous=wrap(densityDiag_lims, alpha=0.1, adjust=2, size=0.25)),
          upper=list(
            continuous=wrap(ggally_cor, stars=F, alpha=1, title="Overall r",
                            digits=2, size=2, justify_labels="left", 
                            fontface="bold", family="mono", colour="grey20")), 
          lower=list(
            continuous=wrap(ggally_points_lims)
          ),
          labeller="label_parsed",
          axisLabels="show") +
  scale_colour_brewer(type="qual", palette="Paired") + 
  scale_fill_brewer(type="qual", palette="Paired") + 
  theme(panel.border=element_rect(colour="grey30", fill=NA)) 
ggsave("figs/pub/Fig_5.png", fig5, width=140, height=130, units="mm", dpi=300)
theme_set(theme_classic())



rank.oos |> 
  filter(.metric %in% c("MCC", "R2-VZ", "PR-AUC", "ROC-AUC", "Schoener's D"),
         covSet %in% c("ensGLM2", "null4wk", paste0("d", 1:16))) |> 
  ungroup() |> 
  left_join(y_i |> select(y, fig_short)) |>
  select(fig_short, model, covSet, PCA, .metric, .estimate) |> 
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  arrange(.metric, fig_short, model) |>
  pivot_wider(names_from=".metric", values_from=".estimate") |> 
  select(-covSet, -PCA, -model) |>
  select(-fig_short) |>
  as.matrix() |>
  cor() |>
  abs() |>
  summary()
  corrplot::corrplot() 

rank.oos |> 
  filter(.metric %in% c("MCC", "R2-VZ", "PR-AUC", "ROC-AUC", "Schoener's D"),
         covSet %in% c("ensGLM2", "null4wk", paste0("d", 1:16))) |> 
  ungroup() |> 
  left_join(y_i |> select(y, fig_short)) |>
  select(fig_short, model, covSet, PCA, .metric, .estimate) |> 
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  arrange(.metric, fig_short, model) |>
  pivot_wider(names_from=".metric", values_from=".estimate") |> 
  select(-covSet, -PCA, -model) |>
  group_by(fig_short) |>
  group_split(.keep=F) |>
  map(~.x |> as.matrix() |> cor()) 
  




# Fig 6 -------------------------------------------------------------------

fig6_mod_labs <- c("Null[Date]", "Ens.", "HBayes", "Ridge", 
                   "MARS", "NN", "RF", "XGB")
fig6_metric_labs <- c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]")

fig6 <- rank.oos  |>
  filter(!grepl("auto|0|Ens-", model)) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=fig6_metric_labs)) |> 
  arrange(y, .metric, rank) |>
  group_by(y, .metric) |>
  mutate(rank=(194-min_rank(rank))/194*100) |>
  ungroup() |>
  ggplot(aes(model, rank, fill=model)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.7, outlier.alpha=0.5, size=0.25) +
  scale_fill_manual("Model", values=mod_cols, guide="none") + 
  labs(x="", y="Rank") +
  scale_x_discrete(labels=parse(text=paste0("italic(", fig6_mod_labs, ")"))) +
  scale_y_continuous("Percentile", breaks=seq(0,100,by=10)) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  theme_ms +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="grey90", linewidth=0.25))

ggsave("figs/pub/Fig_6.png", fig6, width=140, height=90, units="mm", dpi=300)

library(ggdist)
fig6 <- rank.oos  |>
  filter(!grepl("auto|0|Ens-", model)) |>
  # filter(!covSet %in% paste0("d", 13:16)) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=fig6_metric_labs)) |> 
  arrange(y, .metric, rank) |>
  group_by(y, .metric) |>
  mutate(rank=(194-min_rank(rank))/194*100) |>
  ungroup() |>
  ggplot(aes(rank, model, fill=model)) + 
  geom_dots(aes(colour=model), side="bottom", scale=0.5, shape=1,
            binwidth=1, overflow="compress", alpha=0.7) +
  stat_halfeye(scale=0.5, normalize="xy", fill_type="gradient", 
               .width=c(0.66, 0.95), interval_size_range=c(0.2, 0.5),
               fatten_point=0.5,
               aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  scale_colour_manual("Model", values=mod_cols, guide="none") +
  scale_fill_manual("Model", values=mod_cols, guide="none") +
  scale_y_discrete(labels=parse(text=paste0("italic(", fig6_mod_labs, ")"))) +
  scale_x_continuous("Performance percentile", breaks=seq(0,100,by=20)) +
  scale_slab_alpha_continuous(range=c(0.3, 0.9), guide="none") +
  theme_ms + 
  theme(legend.position="none",
        axis.title.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.3),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.15))
ggsave("figs/pub/Fig_6_alt1.png", fig6, width=90, height=120, units="mm", dpi=300)

fig6 <- rank.oos  |>
  filter(!grepl("auto|0|Ens-", model)) |>
  # filter(!covSet %in% paste0("d", 13:16)) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=fig6_metric_labs)) |> 
  arrange(y, .metric, rank) |>
  group_by(y, .metric) |>
  mutate(rank=(194-min_rank(rank))/194*100) |>
  ungroup() |>
  ggplot(aes(rank, .metric, fill=model)) + 
  geom_dots(aes(colour=model), side="bottom", scale=0.5, 
            binwidth=1, overflow="compress", alpha=0.7) +
  stat_halfeye(scale=0.5, normalize="xy", fill_type="gradient", 
               .width=c(0.66, 0.95), interval_size_range=c(0.2, 0.5),
               fatten_point=0.5,
               aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  scale_colour_manual("Model", values=mod_cols, guide="none") +
  scale_fill_manual("Model", values=mod_cols, guide="none") +
  scale_y_discrete(labels=~parse(text=.x), limits=rev(fig6_metric_labs)) +
  scale_x_continuous("Performance percentile", breaks=seq(0,100,by=20)) +
  scale_slab_alpha_continuous(range=c(0.3, 0.9), guide="none") +
  theme_ms + 
  facet_wrap(~model, ncol=1, strip.position="right", labeller=label_parsed) +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.3),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.15))
ggsave("figs/pub/Fig_6_alt2.png", fig6, width=90, height=200, units="mm", dpi=300)

fig6 <- rank.oos  |>
  filter(!grepl("auto|0|Ens-", model)) |>
  # filter(!covSet %in% paste0("d", 13:16)) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=fig6_metric_labs)) |> 
  arrange(y, .metric, rank) |>
  group_by(y, .metric) |>
  mutate(rank=(194-min_rank(rank))/194*100) |>
  ungroup() |>
  left_join(y_i |> select(y, fig_short, fig_long)) |>
  ggplot(aes(rank, fig_short, fill=fig_short)) + 
  geom_dots(aes(colour=fig_short), side="bottom", scale=0.5, 
            binwidth=1, overflow="compress", alpha=0.7) +
  stat_halfeye(scale=0.5, normalize="xy", fill_type="gradient", 
               .width=c(0.66, 0.95), interval_size_range=c(0.2, 0.5),
               fatten_point=0.5,
               aes(slab_alpha=after_stat(-pmax(abs(1-2*cdf), 0.5)))) +
  scale_colour_brewer(type="qual", palette=3, guide="none") +
  scale_fill_brewer(type="qual", palette=3, guide="none") +
  scale_y_discrete() +
  scale_x_continuous("Performance percentile", breaks=seq(0,100,by=20)) +
  scale_slab_alpha_continuous(range=c(0.3, 0.9), guide="none") +
  theme_ms + 
  facet_wrap(~model, ncol=1, strip.position="right", labeller=label_parsed) +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.3),
        panel.grid.minor.x=element_line(colour="grey90", linewidth=0.15))
ggsave("figs/pub/Fig_6_alt3.png", fig6, width=90, height=200, units="mm", dpi=300)



fig6_metrics_df <- rank.oos |>
  filter(!grepl("auto|0|Ens-", model)) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  arrange(y, .metric, rank) |>
  group_by(y, .metric) |>
  mutate(rank=min_rank(rank)) |>
  ungroup()
fig6_df <- fig6_metrics_df |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "MCC", "R2-VZ", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "MCC", "R['VZ']^2", "D[overlap]"))) |> 
  left_join(y_i |> select(y, fig_short, fig_long))
fig6_means <- fig6_df |>
  group_by(model, .metric) |>
  summarise(across(where(is.numeric), mean))
fig6 <- fig6_df |>
  ggplot(aes(model, rank, fill=model)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.7, outlier.alpha=0.5, size=0.25) +
  scale_fill_manual("Model", values=mod_cols, guide="none") + 
  labs(x="", y="Rank") +
  scale_x_discrete(labels=parse(text=paste0("italic(", fig6_mod_labs, ")"))) +
  facet_wrap(~.metric, nrow=1, labeller="label_parsed") + 
  theme_talk +
  theme(legend.position="bottom",
        axis.text.x=element_text(angle=270, hjust=0, vjust=0.5),
        axis.title.x=element_blank(),
        panel.grid=element_blank())

ggsave("figs/pub/Fig_6_TALK.png", fig6, width=160, height=90, units="mm", dpi=300)







# Fig 7 -------------------------------------------------------------------

vi_rng_df <- vi_pt_df |> group_by(fig_short, fig_short_rev, VariableType, varTypeClean) |>
  reframe(mnImp=c(0, max(mnImp)))

fig7 <- vi_pt_df |>
  ggplot(aes(mnImp, fig_short_rev, colour=fig_short)) + 
  geom_line(data=vi_rng_df, linewidth=0.25) +
  geom_point(aes(shape=top, size=top, alpha=top)) + 
  scale_colour_brewer(type="qual", palette=3) + 
  scale_shape_manual(values=c(1, 19), guide="none") +
  scale_size_manual(values=c(0.75, 1.75), guide="none") +
  scale_alpha_manual(values=c(0.3, 1), guide="none") +
  scale_x_continuous("Ensemble weighted relative importance (%)", 
                     breaks=seq(0,100,by=25)) +
  facet_grid(varTypeClean~., scales="free_y", space="free_y", switch="y",
             labeller=label_wrap_gen()) +
  theme_ms +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.9,0.09),
        legend.key.height=unit(0.35, "cm"))
ggsave("figs/pub/Fig_7.png", fig7, width=90, height=200, units="mm", dpi=300)

fig7 <- vi_pt_df |>
  ggplot(aes(fig_short, mnImp, colour=fig_short)) + 
  geom_hline(yintercept=0, colour="grey30", linewidth=0.25) + 
  geom_line(data=vi_rng_df, linewidth=0.25) +
  geom_point(aes(shape=top, size=top, alpha=top)) + 
  scale_colour_brewer(type="qual", palette=3) + 
  scale_shape_manual(values=c(1, 19), guide="none") +
  scale_size_manual(values=c(0.75, 1.75), guide="none") +
  scale_alpha_manual(values=c(0.5, 1), guide="none") +
  scale_y_continuous("Ensemble weighted relative importance (%)", 
                     breaks=seq(0,100,by=25)) +
  facet_grid(.~varTypeClean, scales="free_x", space="free_x", switch="x",
             labeller=label_wrap_gen()) +
  guides(colour=guide_legend(ncol=2)) + 
  theme_talk +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.position="none",
        legend.key.height=unit(0.35, "cm"))
ggsave("figs/pub/Fig_7_TALK.png", fig7, width=260, height=150, units="mm", dpi=300)















# Summaries ---------------------------------------------------------------

# Observed bloom proportion
obs.df |>
  group_by(y) |>
  summarise(propA1=mean(alert=="A1"))
# Bloom starts, continuations, and ends
bloom_df <- obs.df |>
  mutate(year=year(date)) |>
  group_by(y, siteid, year) |>
  mutate(prevAlert=lag(alert),
         bloomNew=(alert=="A1" & (prevAlert=="A0" | is.na(prevAlert))),
         bloomCont=(alert=="A1" & prevAlert=="A1"),
         bloomEnd=(alert=="A0" & prevAlert=="A1")) |>
  group_by(y) |>
  arrange(y, siteid, date) |>
  mutate(bloomID=cumsum(bloomNew)*(alert=="A1"))
# Bloom lengths
bloom_df |> 
  filter(bloomID > 0) |>
  count(bloomID) |>
  summarise(mdLength=median(n),
            mnLength=mean(n),
            q1Length=quantile(n, 0.25),
            q3Length=quantile(n, 0.75),
            q975Length=quantile(n, 0.975),
            varLength=var(n),
            nBlooms=n())
bloom_poisson <- bloom_df |> 
  filter(bloomID > 0) |>
  count(bloomID) |>
  summarise(lambda=mean(n)) |>
  group_by(y) |>
  mutate(pr_bloom=list(map(1:20, ~dpois(.x, lambda)))) |>
  unnest(pr_bloom) |>
  unnest(pr_bloom) |>
  mutate(bloom_length=1:20) |>
  ungroup()

bloom_df |> 
  filter(bloomID > 0) |>
  count(bloomID, name="bloom_length") |>
  group_by(y, bloom_length) |>
  summarise(n_blooms=n()) |>
  group_by(y) |>
  mutate(pr_bloom=n_blooms/sum(n_blooms)) |>
  ungroup() |>
  ggplot(aes(bloom_length, pr_bloom)) + 
  geom_bar(stat="identity") + 
  geom_point(data=bloom_poisson) +
  facet_wrap(~y)
bloom_df |>
  filter(bloomID > 0) |>
  # filter(y=="Alsp") |>
  ggplot(aes(date, siteid, group=bloomID)) + 
  geom_point(shape=1) +
  geom_line(colour="red") +
  facet_wrap(~y, scales="free_y")

# Improvement: ensemble vs. nulls
## IN ORDER OF MANUSCRIPT
### Mean/min/max among all species, metrics
rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(model) |> 
  summarise(gain_q025=quantile(std_gain, 0.025),
            gain_mean=mean(std_gain),
            gain_q975=quantile(std_gain, 0.975))

### Means by monitoring target
rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(y, model) |>
  summarise(across(c("diff_v_null", "std_gain"), 
                   ~paste0(signif(mean(.x),3), " ", signif(min(.x),3), "-", signif(max(.x),3)))) |>
  ungroup() |>
  filter(!grepl("0", model)) |>
  arrange(y, desc(model))
rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(y, model) |>
  summarise(across(c("diff_v_null", "std_gain"), 
                   ~signif(mean(.x),3))) |>
  ungroup() |>
  filter(!grepl("0", model)) |>
  arrange(y, desc(model)) |>
  group_by(y) |>
  summarise(mean_gain=paste0(first(std_gain), " (", last(std_gain), ")"))

### Increase in gain vs. nullDate
rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(y, model) |>
  summarise(across(c("diff_v_null", "std_gain"), mean)) |>
  group_by(y) |>
  filter(!grepl("0", model)) |>
  arrange(desc(model)) |>
  summarise(std_gain_diff=first(std_gain)-last(std_gain),
            std_gain_pctImprove=std_gain_diff/last(std_gain)*100,
            ens_gain=first(std_gain)) |>
  ungroup() |>
  arrange(desc(std_gain_diff)) |>
  summarise(gain_q025=quantile(std_gain_diff, 0.025),
            gain_mean=mean(std_gain_diff),
            gain_q975=quantile(std_gain_diff, 0.975))

### Means by metric
rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=100 * diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(.metric, model) |>
  summarise(across(c("diff_v_null", "std_gain"), 
                   ~paste0(signif(mean(.x),3), " ", signif(min(.x),3), "-", signif(max(.x),3)))) |>
  ungroup() |>
  filter(!grepl("0", model))






# Extra scratch, currently not reported

rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  filter(!grepl("0", model)) |>
  mutate(model=factor(model, levels=c("Null[Date]", "Ensemble"),
                      labels=c("italic(Null[Date])", "italic(Ens.)"))) |>
  ggplot(aes(std_gain, colour=model, fill=model)) + 
  geom_density(alpha=0.5, adjust=0.8) + 
  scale_colour_manual("Model", values=c("grey30", "steelblue3"), 
                      labels=parse_format()) +
  scale_fill_manual("Model", values=c("grey30", "steelblue3"), 
                      labels=parse_format()) +
  labs(x=expression(atop("Standardised gain", 
                         paste("(0: Null" [0]~~"1: Perfect)"))),
       title="Improvement across all monitoring targets and metrics") + 
  xlim(-0.02, 1) +
  theme_talk +
  theme(legend.text.align=0,
        legend.position=c(0.9, 0.7))
ggsave("figs/pub/Fig_2_overview_TALK.png", width=160, height=90, units="mm", dpi=300)

rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  filter(!grepl("0", model)) |>
  arrange(y, .metric, model) |>
  group_by(y, .metric) |>
  summarise(ens_minus_date=last(std_gain)-first(std_gain),
            ens_minus_date_std=(last(std_gain)-first(std_gain))/(1-first(std_gain))) |>
  ggplot(aes(ens_minus_date_std)) + 
  geom_density(aes(colour=y), adjust=1.5, linewidth=0.7) + 
  scale_colour_brewer(type="qual", palette="Paired") +
  geom_density(linewidth=1.25) + 
  xlim(0, 1)

## Mean + range across species within metrics
rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=100 * diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(.metric, model) |>
  summarise(across(c("diff_v_null", "std_gain"), 
                 ~paste0(signif(mean(.x),3), " ", signif(min(.x),3), "-", signif(max(.x),3)))) |>
  ungroup() |>
  filter(!grepl("0", model))

rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(.metric, model) |>
  summarise(across(c("diff_v_null", "std_gain"), mean)) |>
  ungroup() |>
  filter(!grepl("0", model)) |>
  arrange(desc(std_gain))

rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=100 * diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |> 
  filter(!grepl("0", model)) |>
  select(y, model, .metric, std_gain) |>
  pivot_wider(names_from=".metric", values_from="std_gain") |>
  arrange(y, desc(model))



## Mean + range within species across metrics
rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  group_by(y, model) |>
  summarise(across(c("diff_v_null", "std_gain"), mean)) |>
  ungroup() |>
  filter(!grepl("0", model)) |>
  arrange(desc(std_gain))

rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  filter(!grepl("0", model)) |>
  arrange(desc(model), .metric, y) |>
  group_by(y, .metric) |>
  summarise(std_gain_diff=first(std_gain)-last(std_gain),
            std_gain_pctImprove=std_gain_diff/pmax(last(std_gain), 0.001)*100,
            ens_gain=first(std_gain)) |>
  group_by(.metric) |>
  summarise(mean_diff=mean(std_gain_diff))



  


# Candidate model rankings
rank.oos |> 
  filter(!is.na(rank)) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  filter(!grepl("auto|0|Ens-", model)) |>
  arrange(rank) |>
  group_by(y, .metric) |>
  mutate(rank=min_rank(rank)) |>
  mutate(rank=(194-min_rank(rank))/194*100) |>
  group_by(model) |>
  summarise(mdRank=median(rank, type=3),
            q25Rank=quantile(rank, probs=0.25, type=3),
            q75Rank=quantile(rank, probs=0.75, type=3),
            mnRank=mean(rank),
            sdRank=sd(rank),
            seRank=sd(rank)/sqrt(n()),
            iqrRank=IQR(rank),
            MIN=min(rank),
            MAX=max(rank)) |>
  arrange(mnRank) |>
  mutate(across(ends_with("Rank"), ~round(.x, 1)))
  
rank_lm <- aov(rank ~ model * y * .metric, 
              data=rank.oos |> 
                filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
                filter(!grepl("auto|0|Ens-", model)) |>
                arrange(rank) |>
                group_by(y, .metric) |>
                mutate(rank=min_rank(rank)) |>
                ungroup())
TukeyHSD(rank_lm, which="model")$model |>
  as_tibble(rownames="comparison") |>
  filter(grepl("Ensemble|Null", comparison)) |>
  arrange(desc(abs(diff))) 


# Variable importance
vi_pt_df |> 
  mutate(varBase=Variable %>%
           # str_remove_all("Avg") %>%
           str_remove_all("Wk") %>%
           str_remove_all("Delta") %>%
           str_remove_all("Dir[NSEW]") %>%
           str_replace("UXfetchX", "Wind_") %>%
           str_replace("VXfetchX", "Wind_") %>%
           str_replace("UV", "Wind") %>%
           str_replace("U", "Wind") %>%
           str_replace("V", "Wind") %>%
           str_remove_all("1A1") %>%
           str_remove_all("2A1") %>%
           str_remove_all("Wt1") %>%
           str_remove_all("Wt2") %>%
           str_replace("Alert1", "Alert") %>%
           str_replace("Alert2", "Alert") %>%
           str_replace("lnNX", "lnN_") %>%
           str_replace("lonzXlatz", "lonlat") %>%
           str_replace("lonz", "lonlat") %>%
           str_replace("latz", "lonlat") %>%
           str_replace("ydaySinXydayCos", "yday") %>%
           str_replace("ydaySin", "yday") %>%
           str_replace("ydayCos", "yday")) |>
  group_by(varTypeClean, varBase, fig_short) |>
  summarise(totImp=mean(mnImp)) |>
  group_by(varTypeClean, varBase) |>
  mutate(aggTot=sum(totImp)) |>
  ungroup() |>
  arrange(aggTot) |>
  mutate(varBase_agg=factor(varBase, levels=unique(varBase))) |>
  ggplot(aes(totImp, varBase_agg, fill=fig_short)) + 
  geom_bar(stat="identity") + 
  scale_fill_brewer(type="qual", palette=3) +
  facet_grid(varTypeClean~., scales="free", space="free_y", switch="y",
             labeller=label_wrap_gen())








# Fig S1 ------------------------------------------------------------------

figS1 <- vi_pt_df |> 
  filter(top) |> 
  mutate(Variable=str_remove(Variable, "Dir[NSEW]") |>
           str_remove_all("Wk") |>
           str_remove_all("Wt") |>
           str_replace_all("UV", "Wind") |>
           str_replace_all("UX", "WindX") |>
           str_replace_all("VX", "WindX")) |>
  mutate(varTypeClean=forcats::fct_recode(varTypeClean, "Prev. year"="Previous year")) |>
  ggplot(aes(y, Variable, colour=y)) + 
  geom_point() + 
  scale_colour_brewer(type="qual", palette="Paired") +
  facet_grid(varTypeClean~., scales="free_y", space="free_y",
             labeller=label_wrap_gen(8)) +
  theme_ms +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.key.height=unit(0.35, "cm"),
        panel.grid.major.y=element_line(colour="grey", linewidth=0.1))
ggsave("figs/pub/Fig_S1.png", figS1, width=90, height=240, units="mm", dpi=300)







# Fig S2 ------------------------------------------------------------------

figS2 <- rank.oos |>
  filter(covSet %in% c("nullGrand", "null4wk", "ensGLM2")) |>
  filter(.metric %in% c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D")) |>
  group_by(.metric, y) |>
  arrange(model) |>
  mutate(diff_v_null=.estimate - first(.estimate),
         std_gain=diff_v_null / 
           if_else(.metric=="Schoener's D", -first(.estimate), (1-first(.estimate)))) |>
  mutate(diff_v_null=diff_v_null * if_else(.metric=="Schoener's D", -1, 1)) |>
  filter(!grepl("0", model)) |>
  select(y, model, .metric, std_gain) |>
  left_join(y_i |> select(y, fig_short)) |>
  pivot_wider(names_from="model", values_from="std_gain") |>
  mutate(.metric=factor(.metric, 
                        levels=c("ROC-AUC", "PR-AUC", "R2-VZ", "MCC", "Schoener's D"),
                        labels=c("ROC-AUC", "PR-AUC", "R['VZ']^2", "MCC", "D[overlap]"))) |>
  ggplot(aes(`Null[Date]`, Ensemble, colour=fig_short)) + 
  geom_abline(linewidth=0.25) + 
  geom_point(aes(shape=.metric), size=2) + 
  geom_rug(linewidth=0.2, length=unit(0.2, "cm")) +
  scale_colour_brewer("Monitoring target", type="qual", palette=3, 
                      guide=guide_legend(ncol=2, order=1)) + 
  scale_shape_manual("Performance metric", values=c(1, 19, 3, 2, 4), 
                     labels=parse_format()) +
  scale_x_continuous(expression(italic(Null['Date'])~'standardised gain'), 
                     breaks=c(0, 0.5, 1), limits=c(-0.05,1)) + 
  scale_y_continuous(expression(italic(Ensemble)~"standardised gain"), 
                     breaks=c(0, 0.5, 1), limits=c(-0.05,1)) + 
  theme_ms +
  theme(legend.position=c(0.855, 0.325),
        # legend.spacing.y=unit(0.05, "cm"),
        legend.margin=margin(-1,0,0,0),
        legend.text.align=0,
        legend.key.height=unit(0.35, "cm"),
        legend.key.width=unit(0.35, "cm"))
ggsave("figs/pub/Fig_S2.png", figS2, width=90, height=85, units="mm", dpi=300)








# Fig S4 ------------------------------------------------------------------

# Spatial metrics by site
out.oos |>
  filter(covSet == "ensGLM2") |>
  select(y, siteid, model, alert, predMCC) |>
  group_by(y, model, siteid) |>
  na.omit() |>
  mcc(predMCC, truth=alert) |>
  filter(!is.na(.estimate)) |>
  left_join(y_i |> select(y, fig_long, type)) |>
  left_join(bind_rows(site_tox.df |> select(siteid, lon, lat) |> mutate(type="tox"),
                      site_hab.df |> select(siteid, lon, lat) |> mutate(type="hab"))) |>
  ggplot() + 
  geom_sf(data=coast) +
  geom_point(aes(lon, lat, colour=.estimate), size=2) + 
  scale_colour_gradient2("MCC") +
  facet_wrap(~fig_long, nrow=2) + 
  theme(legend.position=c(0.9, 0.2))




out.oos |>
  filter(covSet %in% c("ensGLM2", "null4wk", "nullGrand")) |>
  mutate(month=month(date, label=TRUE)) |>
  select(y, model, covSet, PCA, month, alert, prA1) |>
  na.omit() |>
  mutate(nA1=sum(alert=="A1")) |>
  ungroup() |>
  filter(nA1 > 2) |>
  group_by(y, month, model) |>
  find_AUCPR_min(y) |>
  nest(dat=c(prA1, alert)) |>
  mutate(AUCPR=map_dbl(dat, ~average_precision(.x, alert, prA1, event_level="second")$.estimate),
         AUCNPR=(AUCPR-AUCPR_min)/(1-AUCPR_min)) |>
  select(-dat) |>
  rename(.estimate=AUCNPR) |>
  ggplot(aes(month, .estimate, colour=model, group=model)) + 
  geom_line() + 
  facet_wrap(~y) + 
  ylab("PR-AUC")

out.oos |>
  filter(covSet %in% c("ensGLM2", "null4wk", "nullGrand")) |>
  mutate(month=month(date, label=TRUE)) |>
  select(y, month, model, alert, prA1) |>
  group_by(y, model, month) |>
  na.omit() |>
  mutate(nA1=sum(alert=="A1")) |>
  ungroup() |>
  filter(nA1 > 2) |>
  group_by(y, month, model) |>
  roc_auc(prA1, truth=alert, event_level="second") |>
  filter(!is.na(.estimate)) |>
  ggplot(aes(month, .estimate, colour=model, group=model)) + 
  geom_line() + 
  facet_wrap(~y) + 
  ylab("ROC-AUC")

out.oos |>
  filter(covSet %in% c("ensGLM2", "null4wk", "nullGrand")) |>
  mutate(month=month(date, label=TRUE)) |>
  select(y, month, model, alert, predMCC) |>
  group_by(y, model, month) |>
  na.omit() |>
  mutate(nA1=sum(alert=="A1")) |>
  ungroup() |>
  filter(nA1 > 2) |>
  group_by(y, month, model) |>
  mcc(predMCC, truth=alert) |>
  filter(!is.na(.estimate)) |>
  ggplot(aes(month, .estimate, colour=model, group=model)) + 
  geom_hline(yintercept=0, linetype=3) + 
  geom_line() + 
  facet_wrap(~y) + 
  ylab("MCC")


out.oos |>
  filter(covSet %in% c("ensGLM2", "null4wk", "nullGrand")) |>
  mutate(month=month(date, label=TRUE)) |>
  mutate(prA1=if_else(prA1==0, 1e-5, prA1),
         prA1=if_else(prA1==1, 1-1e-5, prA1),
         alert=as.numeric(alert=="A1")) |>
  group_by(y, month, model, PCA, covSet) |>
  mutate(nA1=sum(alert==1)) |>
  ungroup() |>
  filter(nA1 > 2) |>
  calc_R2(type="vz", y, month) |>
  rename(.estimate=R2) |>
  mutate(.estimate=pmin(pmax(.estimate, 0), 1)) |>
  ggplot(aes(month, .estimate, colour=model, group=model)) + 
  geom_line() + 
  facet_wrap(~y) + 
  ylab("R2-VZ")

library(kerneval)
out.oos |>
  filter(covSet %in% c("ensGLM2", "null4wk")) |>
  mutate(month=month(date, label=TRUE)) |>
  select(month, y, covSet, PCA, model, obsid, alert, prA1) %>%
  filter(!is.na(prA1)) |>
  group_by(y, month, model, PCA, covSet) |>
  mutate(nA1=sum(alert=="A1")) |>
  ungroup() |>
  filter(nA1 > 2) |>
  pivot_wider(names_from="alert", values_from="prA1") |>
  group_by(y, month, model, PCA, covSet) |>
  summarise(.estimate=schoenr(density(A0, na.rm=T), density(A1, na.rm=T))) |>
  ggplot(aes(month, .estimate, colour=y, group=paste(y, model), shape=model)) + 
  geom_point() +
  scale_shape_manual(values=c(1, 19)) +
  scale_colour_brewer(type="qual", palette="Paired") +
  ylab("D")







# Null date ---------------------------------------------------------------

out.oos |> 
  filter(model %in% c("Null[Date]", "Null[0]")) |> 
  mutate(yday=yday(date)) |> 
  group_by(yday, y, model) |> 
  summarise(prA1=mean(prA1)) |> 
  ungroup() |>
  mutate(date=ymd("2020-12-31") + yday,
         panel=case_when(y %in% c("Alsp", "PSP") ~ "Alexandrium, PST",
                         y %in% c("Disp", "DSP") ~ "Dinophysis, DST",
                         y %in% c("Pssp", "Psde", "Psse", "ASP") ~ "Pseudo-nitzschia, AST",
                         y %in% c("Kami") ~ "Karenia mikimotoi")) |>
  left_join(y_i |> select(y, fig_short)) |>
  ggplot(aes(date, prA1, colour=fig_short, linetype=model)) + 
  geom_vline(xintercept=ymd(paste0("2021-", c("01", "04", "07", "10"), "-01")),
             linewidth=0.25, colour="grey90") +
  geom_hline(yintercept=0, linewidth=0.25, colour="grey90") +
  geom_line(linewidth=0.7) + 
  scale_colour_brewer(type="qual", palette="Paired") + 
  scale_linetype_manual(values=c(2, 1)) + 
  scale_x_date("", date_breaks="1 month", date_labels="%b") + 
  ylab("Null model alert probabilities") + 
  facet_wrap(~panel) + 
  theme_talk + 
  theme(legend.position="none", 
        axis.title.x=element_blank())
ggsave("figs/pub/Fig_0_NullDate_TALK.png", width=200, height=120, units="mm", dpi=300)







# Shellfish sites ---------------------------------------------------------

mss_sites <- read_csv("data/ms_site_details.csv") |>
  janitor::clean_names("small_camel") |>
  filter(aquacultureType=="Shellfish") |>
  filter(producingInLast3Years=="Yes")
write_csv(mss_sites, "data/mss_sites_shellfish.csv")



# PR curves ---------------------------------------------------------------

thresholds <- out.oos |> 
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  select(y, obsid, model, prevAlert, optMCC) |> 
  group_by(model, y, prevAlert) |>
  summarise(thresh_MCC=mean(optMCC))

out.oos |> 
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  # na.omit() |>
  group_by(y, model, PCA, covSet) |> 
  pr_curve(prA1, truth=alert, event_level="second") |> 
  filter(recall > 0) |>
  ggplot(aes(recall, precision, colour=model, group=paste(PCA, covSet, model))) + 
  geom_line() + 
  scale_colour_manual(values=c("grey", "black"), labels=parse_format()) +
  facet_wrap(~y) + ggtitle("oos") + 
  theme(legend.text.align=0)
out.oos |> 
  filter(covSet %in% c("null4wk", "ensGLM2")) |>
  # na.omit() |>
  group_by(y, model, PCA, covSet) |> 
  roc_curve(prA1, truth=alert, event_level="second") |> 
  ggplot(aes(1-specificity, sensitivity, colour=model, group=paste(PCA, covSet, model))) + 
  geom_abline(linetype=3, colour="grey") + 
  geom_line() + 
  scale_colour_manual(values=c("grey", "black"), labels=parse_format()) +
  facet_wrap(~y) + ggtitle("oos") + 
  theme(legend.text.align=0)









# Anim S1 -----------------------------------------------------------------


site_month_sum.df <- site_all.df |> select(-type) |>
  inner_join(obs.mo.Site) |>
  mutate(month=month(dateStd, label=T, abbr=F))

for(i in 1:12) {
  p <- site_month_sum.df |> 
    filter(month==month.name[i]) |>
    ggplot() + 
    geom_sf(data=coast.sf, fill="grey50", colour="grey50") + 
    geom_point(aes(lon, lat, colour=prA, group=sin), shape=1) + 
    scale_x_continuous("Longitude", breaks=c(-8, -4)) +
    scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
    scale_colour_viridis_c("Proportion of\nrecords ≥ amber\n", option="plasma", limits=c(0,1), end=0.95) +
    # scale_colour_viridis_c("(warnings+alerts)/N\n", option="plasma", limits=c(0,1), end=0.95) +
    ggtitle(month.name[i]) +
    theme_ms +
    theme(legend.position=c(0.9, 0.2),
          legend.key.width=unit(1, "mm"),
          panel.background=element_rect(colour=NA, fill="grey70")) + 
    facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
  ggsave(glue("figs/temp/Anim_S1_{str_pad(i, 2, 'left', '0')}.png"), p, 
            width=190, height=145, dpi=300, units="mm")
}

av::av_encode_video(dirf("figs/temp", "Anim_S1_"), 
                    glue("figs/pub/Anim_S1.mp4"),
                    framerate=1)

for(i in 1:12) {
  p <- site_month_sum.df |> 
    filter(month==month.name[i]) |>
    arrange(prA) |>
    ggplot() + 
    geom_sf(data=coast.sf, fill="grey50", colour="grey50") + 
    geom_point(aes(lon, lat, colour=prA*100, group=sin), shape=1) + 
    scale_x_continuous("Longitude", breaks=c(-6, -2)) +
    scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
    scale_colour_viridis_c("% Warn\nor Alert\n", option="plasma", limits=c(0,100), end=0.95) +
    ggtitle(month.name[i]) +
    theme_talk +
    theme(legend.position=c(0.9, 0.22),
          legend.key.width=unit(2, "mm"),
          axis.title=element_blank(),
          panel.background=element_rect(colour=NA, fill="grey70")) + 
    facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
  ggsave(glue("figs/temp/Anim_TALK_S1_{str_pad(i, 2, 'left', '0')}.png"), p, 
         width=190, height=145, dpi=300, units="mm")
}

av::av_encode_video(dirf("figs/temp", "Anim_TALK_S1_"), 
                    glue("figs/pub/Anim_S1_TALK.mov"),
                    framerate=1)
sevcheck::img_to_gif(dirf("figs/temp", "Anim_TALK_S1_"),
                     2,
                     glue("figs/pub/Anim_S1_TALK.gif"))



grd_df <- readRDS("data/grid_scot_coast.rds") |>
  filter(date==first(date)) |>
  st_as_sf(coords=c("lon", "lat"), crs=27700, remove=F)
site_grd_df <- site_all.df |>
  st_as_sf(coords=c("lon", "lat"), crs=27700) %>%
  mutate(grd_i=st_nearest_feature(., grd_df))
nbr_month_sum.df <- site_month_sum.df |>
  left_join(site_grd_df %>% st_drop_geometry() |> select(sin, grd_i), by="sin") |>
  group_by(y, fig_long, month, grd_i) |>
  summarise(prA=sum(prA*N)/sum(N)) |>
  ungroup() |>
  mutate(lon=grd_df$lon[grd_i],
         lat=grd_df$lat[grd_i])

for(i in 1:12) {
  p <- nbr_month_sum.df |>
    filter(month==month.name[i]) |>
    arrange(prA) |>
    ggplot() + 
    geom_sf(data=coast.sf, fill="grey30", colour="grey30") + 
    geom_point(aes(lon, lat, colour=prA*100), shape=1) + 
    scale_x_continuous("Longitude", breaks=c(-6, -2)) +
    scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
    scale_colour_viridis_c("% Warn\nor alert\n", option="plasma", limits=c(0,100), end=0.95) +
    # ggtitle(month.name[i]) +
    annotate("text", x=20000, y=1200000, label=month.abb[i], colour="grey90", hjust=0) +
    theme_talk +
    theme(legend.position=c(0.9, 0.22),
          legend.key.width=unit(2, "mm"),
          axis.title=element_blank(),
          panel.background=element_rect(colour=NA, fill="grey50")) + 
    facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
  ggsave(glue("figs/temp/Anim_TALK_S1_{str_pad(i, 2, 'left', '0')}.png"), p, 
         width=210, height=145, dpi=300, units="mm")
}

av::av_encode_video(dirf("figs/temp", "Anim_TALK_S1_"), 
                    glue("figs/pub/Anim_S1_TALK.mov"),
                    framerate=1)
sevcheck::img_to_gif(dirf("figs/temp", "Anim_TALK_S1_"),
                     2,
                     glue("figs/pub/Anim_S1_TALK.gif"))






# Anim S2 -----------------------------------------------------------------

# coast_lines.sf <- coast.sf |>
#   st_cast("MULTILINESTRING")
# coast_innerBuffer.sf <- coast.sf %>%
#   st_intersection(st_buffer(coast_lines.sf, 5e3), .) |>
#   st_sym_difference(coast.sf)
# coast_buffer <- coast_lines.sf %>%
#   st_buffer(25e3) %>%
#   st_union()
# 
# dat <- readRDS("data/0_init/ASP_16-Avg1_Xf1_XN1_Del1_dataSplit.rds") |>
#   training() 
# grd_df <- expand_grid(date=ymd("2021-01-01") + seq(0, 362, by=7),
#                       lon=seq(min(dat$lon)-50e3, max(dat$lon)+50e3, by=5e3), 
#                       lat=seq(min(dat$lat)-50e3, max(dat$lat)+50e3, by=5e3)) |>
#   group_by(date) |> 
#   mutate(id=row_number()) |> ungroup() |>
#   mutate(yday=yday(date))
# grd_sf <- grd_df |> 
#   filter(yday==1) |>
#   st_as_sf(coords=c("lon", "lat"), remove=F, crs=27700) |>
#   st_intersection(coast_buffer) |>
#   st_difference(coast_innerBuffer.sf) |>
#   filter(lat > min(lat))
# grd_df |> 
#   filter(id %in% grd_sf$id) |> 
#   rename(siteid=id) |>
#   mutate(obsid=row_number()) |>
#   saveRDS("data/grid_scot_coast.rds")

spatTime_df <- readRDS("out/clean/out_spatTime.rds") |>
  filter(model=="Ensemble2", yday != max(yday)) |>
  mutate(date=ymd("2020-12-31") + yday) |>
  ungroup() |>
  left_join(y_i |> select(abbr, type, fig_long), by=c("y"="abbr")) |>
  select(y, date, lon, lat, prevAlert, alert, yday, model, prA1, fig_long)

library(gganimate)
wk.def <- tibble(yday=1:364,
                 wkNum=rep(1:26, each=14))
obs.wk <- obs.df |>
  mutate(yday=yday(date)) |>
  left_join(wk.def) |>
  group_by(y, wkNum) |>
  summarise(pAlertAvg=mean(alert=="A1")) |>
  ungroup()

S2_df <- spatTime_df |>
  left_join(wk.def) |>
  left_join(obs.wk) |>
  mutate(prA1_wt=prA1 * if_else(prevAlert=="A0", 1-pAlertAvg, pAlertAvg)) |>
  group_by(y, wkNum, date, lon, lat, model, fig_long) |>
  summarise(prA1=sum(prA1_wt)) |>
  ungroup()
pr_rng <- range(S2_df$prA1, na.rm=T)

for(i in 1:(n_distinct(spatTime_df$date)-1)) {
  i_df <- S2_df |>
    filter(wkNum == i) 
  p <- i_df |>
    ggplot() + 
    geom_tile(aes(lon, lat, fill=prA1)) +
    scale_fill_viridis_c("Spatiotemporal\nconditional\nprobability", limits=c(0, pr_rng[2]), option="plasma") +
    geom_sf(data=coast.sf, colour=NA) +
    ggtitle(glue("Weighted mean state[t-1]: {month(first(i_df$date), label=T)}-{day(first(i_df$date))}")) +
    theme_ms +
    theme(legend.position=c(0.9, 0.2),
          legend.key.width=unit(1, "mm"),
          axis.text=element_blank(),
          axis.title=element_blank()) + 
    facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
  ggsave(glue("figs/temp/anim_S2_{str_pad(i, 2, 'left', '0')}.png"), 
         p, width=190, height=150, dpi=400, units="mm")
}
av::av_encode_video(dirf("figs/temp", "anim_S2"), 
                    glue("figs/pub/Anim_S2.mp4"),
                    framerate=1)





# Anim S3 -----------------------------------------------------------------

max_pr <- spatTime_df |>
  left_join(wk.def) |>
  left_join(obs.wk) |>
  mutate(prA1_wt=prA1 * if_else(prevAlert=="A0", 1-pAlertAvg, pAlertAvg)) |>
  group_by(y, wkNum, date, lon, lat, model, fig_long) |>
  summarise(prA1=sum(prA1_wt)) |>
  ungroup() |>
  filter(wkNum < 26) |>  # no observations for end Dec
  group_by(y, model) |>
  mutate(prA1_rel=(prA1-min(prA1))/(max(prA1)-min(prA1))) |>
  ungroup() |>
  group_by(fig_long, model) |>
  summarise(lon=min(lon), 
            lat=max(lat), 
            pr=paste0("Max pr(alert):\n", round(max(prA1), 3)))
S3_df <- spatTime_df |>
  left_join(wk.def) |>
  left_join(obs.wk) |>
  mutate(prA1_wt=prA1 * if_else(prevAlert=="A0", 1-pAlertAvg, pAlertAvg)) |>
  group_by(y, wkNum, date, lon, lat, model, fig_long) |>
  summarise(prA1=sum(prA1_wt)) |>
  ungroup() |>
  filter(wkNum < 26) |>  # no observations for end Dec
  group_by(y, model) |>
  mutate(prA1_rel=(prA1-min(prA1))/(max(prA1)-min(prA1))) |>
  ungroup()

for(i in 1:n_distinct(S3_df$date)) {
  i_df <- S3_df |>
    filter(wkNum == i) 
  p <- i_df |>
    ggplot() + 
    geom_tile(aes(lon, lat, fill=prA1_rel)) +
    scale_fill_viridis_c("Spatiotemporal\nrelative\nconditional\nprobability",
                         limits=c(0, 1), option="plasma") +
    geom_sf(data=coast.sf, colour=NA) +
    geom_text(data=max_pr, x=7011, y=1235757, aes(label=pr), hjust=0, vjust=1, size=2) +
    ggtitle(glue("Weighted mean state[t-1]: {month(first(i_df$date), label=T)}-{day(first(i_df$date))}")) +
    theme_ms +
    theme(legend.position=c(0.9, 0.2),
          legend.key.width=unit(1, "mm"),
          axis.text=element_blank(),
          axis.title=element_blank()) + 
    facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
  ggsave(glue("figs/temp/anim_S3_{str_pad(i, 2, 'left', '0')}.png"), 
         p, width=190, height=150, dpi=400, units="mm")
}
av::av_encode_video(dirf("figs/temp", "anim_S3"), 
                    glue("figs/pub/Anim_S3.mp4"),
                    framerate=1)









# Anim S4 -----------------------------------------------------------------

obs_week <- expand_grid(date_week=paste0(c(2015:2022), "-01-01") |> 
                          ymd() |> map(~.x + (seq(1,52,by=2))*7) |> list_c(),
                        siteid=1:max(obs.df$siteid),
                        y=unique(obs.df$y)) |>
  left_join(obs.df |>
              mutate(year=year(date),
                     week=floor(week(date)/2)*2,
                     date_week=ymd(paste0(year, "-01-01")) + (week-1)*7) |>
              select(y, siteid, date_week, alert)) |>
  group_by(y, siteid, date_week) |>
  summarise(alert=if_else(any(alert=="A1"), "1: Yes", "0: No")) |>
  ungroup() |>
  left_join(y_i |> select(abbr, fig_long, type) |> rename(y=abbr)) |>
  right_join(bind_rows(site_hab.df |> select(siteid, sin, lon, lat) |> mutate(type="hab"),
                       site_tox.df |> select(siteid, sin, lon, lat) |> mutate(type="tox"))) |>
  filter(!is.na(alert))
animS4 <- ggplot(coast.sf) + 
  geom_sf(fill="grey80", colour="grey80") + 
  geom_point(data=obs_week, aes(lon, lat, colour=alert, shape=alert, group=sin),) + 
  transition_manual(date_week) +
  scale_x_continuous("Longitude", breaks=c(-8, -4)) +
  scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
  scale_colour_manual("Warning or alert", values=c("grey40", "red2"), na.value=NA) +
  scale_shape_manual("Warning or alert", values=c(1, 19), na.value=NA) +
  ggtitle("{current_frame}") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm"),
        axis.title=element_blank(),
        axis.text=element_blank()) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
anim_save("figs/pub/anim_S4_2fps.mp4", animS4, nframes=n_distinct(obs_week$date_week), 
          fps=2, renderer=av_renderer(), 
          width=190, height=145, res=400, units="mm")








# Anim S5 -----------------------------------------------------------------

pred_week <- expand_grid(date_week=paste0(c(2016:2022), "-01-01") |> 
                           ymd() |> map(~.x + (seq(1,52,by=2))*7) |> list_c(),
                         siteid=1:max(obs.df$siteid),
                         y=unique(obs.df$y)) |>
  left_join(out.oos |>
              filter(covSet=="ensGLM2") |>
              mutate(year=year(date),
                     week=floor(week(date)/2)*2,
                     date_week=ymd(paste0(year, "-01-01")) + (week-1)*7) |>
              select(y, siteid, date_week, alert, prA1)) |>
  group_by(y, siteid, date_week) |>
  summarise(alert=if_else(any(alert=="A1"), "A1", "A0"),
            prA1=mean(prA1, na.rm=T)) |>
  ungroup() |>
  left_join(y_i |> select(abbr, fig_long, type) |> rename(y=abbr)) |>
  right_join(bind_rows(site_hab.df |> select(siteid, sin, lon, lat) |> mutate(type="hab"),
                       site_tox.df |> select(siteid, sin, lon, lat) |> mutate(type="tox"))) |>
  filter(!is.na(alert))
animS5 <- ggplot(coast.sf) + 
  geom_sf(fill="grey80", colour="grey80") + 
  geom_point(data=pred_week |> filter(year(date_week)==2018), aes(lon, lat, colour=prA1, group=sin),) + 
  transition_manual(date_week) +
  scale_x_continuous("Longitude", breaks=c(-8, -4)) +
  scale_y_continuous("Latitude", breaks=c(55, 58, 61)) +
  scale_colour_viridis_c("Forecasted P(alert)", option="turbo", na.value=NA) +
  ggtitle("{current_frame}") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm"),
        axis.title=element_blank(),
        axis.text=element_blank()) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
anim_save("figs/pub/anim_S5_2fps.mp4", animS5, nframes=25,#n_distinct(pred_week$date_week), 
          fps=2, renderer=av_renderer(), 
          width=190, height=145, res=400, units="mm")







# deprecated --------------------------------------------------------------


animS2a <- spatTime_df |>
  filter(prevAlert=="A0") |>
  ggplot() + 
  geom_tile(aes(lon, lat, fill=prA1)) +
  scale_fill_viridis_c("Spatiotemporal\nconditional\nprobability", limits=c(0, NA), option="inferno") +
  geom_sf(data=coast.sf, colour=NA) +
  transition_manual(date) + 
  ggtitle("t-1 = No alert: {month(current_frame, label=T)}-{day(current_frame)}") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm"),
        axis.text=element_blank()) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
anim_save("figs/pub/anim_S2a.mp4", animS2a, nframes=n_distinct(spatTime_df$date), fps=5, 
          width=190, height=145, res=400, units="mm")
animS2b <- spatTime_df |>
  filter(prevAlert=="A1") |>
  ggplot() + 
  geom_tile(aes(lon, lat, fill=prA1)) +
  scale_fill_viridis_c("Spatiotemporal\nconditional\nprobability", limits=c(0, NA), option="inferno") +
  geom_sf(data=coast.sf, colour=NA) +
  transition_manual(date) + 
  ggtitle("t-1 = Alert: {month(current_frame, label=T)}-{day(current_frame)}") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm"),
        axis.text=element_blank()) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
anim_save("figs/pub/anim_S2b.mp4", animS2b, nframes=n_distinct(spatTime_df$date), fps=5, 
          width=190, height=145, res=400, units="mm")
animS3a <- spatTime_df |>
  filter(prevAlert=="A0") |>
  ggplot() + 
  geom_tile(aes(lon, lat, fill=prA1_rel)) +
  scale_fill_viridis_c("Spatiotemporal\nrelative\nconditional\nprobability", limits=c(0, 1), option="inferno") +
  geom_sf(data=coast.sf, colour=NA) +
  transition_manual(date) +
  ggtitle("t-1 = No alert: {month(current_frame, label=T)}-{day(current_frame)}") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm"),
        axis.text=element_blank()) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
anim_save("figs/pub/anim_S3a.mp4", animS3a, nframes=n_distinct(spatTime_df$date), 
          fps=5, renderer=av_renderer(), 
          width=190, height=145, res=400, units="mm")
animS3b <- spatTime_df |>
  filter(prevAlert=="A1") |>
  ggplot() + 
  geom_tile(aes(lon, lat, fill=prA1_rel)) +
  scale_fill_viridis_c("Spatiotemporal\nrelative\nconditional\nprobability", limits=c(0, 1), option="inferno") +
  geom_sf(data=coast.sf, colour=NA) +
  transition_manual(date) +
  ggtitle("t-1 = Alert: {month(current_frame, label=T)}-{day(current_frame)}") +
  theme_ms +
  theme(legend.position=c(0.9, 0.2),
        legend.key.width=unit(1, "mm"),
        axis.text=element_blank()) + 
  facet_wrap(~fig_long, nrow=2, labeller=label_wrap_gen(17))
anim_save("figs/pub/anim_S3b.mp4", animS3b, nframes=n_distinct(spatTime_df$date), 
          fps=5, renderer=av_renderer(), 
          width=190, height=145, res=400, units="mm")











# taylor diagram ----------------------------------------------------------


mods <- c("XGB", "HB", "MARS", "NN", "RF", "Ridge", "Null[Date]", "Ensemble")
tay_cols <- mod_cols
tay_cols[1:11] <- "black"
par(mfrow=c(2,5))
for(y_t in levels(out.oos$y)) {
  oos_y <- out.oos |> 
    filter(y==y_t) |>
    mutate(alertNum=as.numeric(alert))
  oos_null <- oos_y |>
    filter(model=="Null[0]")
  taylor.diagram(oos_null$alertNum, 
                 oos_null$prA1, 
                 pch=4,
                 col=tay_cols[1], 
                 main=y_i$fig_long[y_i$abbr==y_t],
                 ref.sd=T)
  
  for(m in mods) {
    oos_m <- oos_y |>
      filter(model==m)
    runs <- unique(oos_m$run)
    for(r in runs) {
      taylor.diagram(filter(oos_m, run==r)$alertNum, 
                     filter(oos_m, run==r)$prA1,
                     col=tay_cols[[m]], 
                     add=T, 
                     pcex=ifelse(grepl("Null|Ensemble", m), 2, 0.75),
                     pch=switch(m,
                                "Null[Date]"=8,
                                "Ensemble"=19,
                                "XGB"=,
                                "HB"=,
                                "MARS"=,
                                "NN"=,
                                "RF"=,
                                "Ridge"=1))
    }
  }
  
  
}


