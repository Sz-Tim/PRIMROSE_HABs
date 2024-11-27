
library(tidyverse)
library(sevcheck)
library(glue)

out.dir <- "out/0_init"

y_i <- bind_rows(read_csv("data/i_hab.csv") |> arrange(abbr) |> mutate(type="hab"),
                 read_csv("data/i_tox.csv") |> arrange(abbr) |> mutate(type="tox")) |>
  mutate(fig_short=factor(fig_short, 
                          levels=c("A. sp.", "PSTs", "D. sp.", "DSTs", 
                                   "P. sp.", "P. del.", "P. ser.", "DA", 
                                   "K. mik.")),
         fig_long=factor(fig_long, 
                         levels=c("Alexandrium", "PSTs", "Dinophysis", "DSTs (OA/DTXs/PTXs)",
                                  "Karenia mikimotoi",
                                  "Pseudo-nitzschia", 
                                  "Pseudo-nitzschia delicatissima", "Pseudo-nitzschia seriata", 
                                  "Domoic acid")),
         y=factor(abbr, levels=c("Alsp", "PSP", "Disp", "DSP", "Pssp", 
                                 "Psde", "Psse", "ASP", "Kami"))) |>
  drop_na()
col_cmems <- readRDS("data/cmems_vars.rds")
col_wrf <- readRDS("data/wrf_vars.rds")
varTypes <- list(
  spacetime=c("yday", "ydayCos", "ydaySin", "ydaySinXydayCos",
              "lon", "lat", "latz", "lonz", "lonzXlatz"),
  autoreg=c("lnNWt1", "lnNAvg1", "prAlertAvg1", "alert1A1",
            "lnNWt2", "lnNAvg2", "prAlertAvg2", "alert2A1"),
  prevYr=c("lnNPrevYr", "lnNAvgPrevYr", "prAlertPrevYr", "prAlertAvgPrevYr"),
  fetch="fetch",
  CMEMS=col_cmems, 
  WRF=col_wrf,
  UV_interact=c(
    paste("UWkXfetch", grep("Dir", col_cmems, value=T), sep="X"),
    paste("VWkXfetch", grep("Dir", col_cmems, value=T), sep="X"),
    paste("UWkXfetch", grep("^[Precip|Shortwave|sst].*Dir", col_wrf, value=T), sep="X"),
    paste("VWkXfetch", grep("^[Precip|Shortwave|sst].*Dir", col_wrf, value=T), sep="X")
  ),
  hab=c(outer(filter(y_i, type=="hab")$abbr, c("lnNAvg", "prA"), "paste0"))
)
varTypes$auto_interact <- paste("lnNWt1", 
                                c(varTypes$autoreg[-1], varTypes$prevYr, 
                                  varTypes$fetch, varTypes$CMEMS, varTypes$WRF), 
                                sep="X")

varType_df <- imap_dfr(varTypes, ~tibble(VariableType=.y, variable=.x)) |>
  mutate(varTypeClean=factor(VariableType, 
                             levels=names(varTypes),
                             labels=c("Space time", "Autoregression",
                                      "Previous year", "Fetch", "Biogeochemistry",
                                      "Weather", "Wind interaction",
                                      "HAB", "Autoregression interaction"))) |>
  mutate(varTypeClean=lvls_reorder(varTypeClean, c(1,2,9,3:8))) |>
  arrange(varTypeClean, variable) |>
  mutate(variable_ordered=factor(variable, levels=unique(variable)))

for(i in 1:nrow(y_i)) {
  y <- y_i$y[i]
  try({
    
  adp_df <- dirrf(glue("{out.dir}/compiled"), glue("{y}.*_adp")) |> 
    map(readRDS) |> 
    reduce(full_join, by=join_by(variable, x)) |> 
    pivot_longer(-(1:2)) |>
    drop_na() |>
    mutate(model=str_split_fixed(
      str_remove(name, "PCA.") |> str_remove("_alert_A1"), "\\.", 2)[,2]) |>
    inner_join(varType_df, by=join_by(variable))
  adp_ens <- readRDS(glue("out/0_init/compiled/{y}_adp_Ensemble.rds")) |>
    inner_join(varType_df) |>
    rename(value=yhat)
  
  p <- adp_df |>
    ggplot(aes(x, value, colour=varTypeClean)) + 
    geom_hline(yintercept=0, linewidth=0.4, colour="black") +
    geom_line(aes(group=name), alpha=0.5, linewidth=0.2) + 
    geom_line(data=adp_ens, linewidth=0.5, colour="black") +
    scale_colour_brewer(type="qual", palette=2, guide="none") +
    {if(y_i$type[i]=="hab") facet_wrap(~variable_ordered, scales="free_x", nrow=14)} +
    {if(y_i$type[i]=="tox") facet_wrap(~variable_ordered, scales="free_x", nrow=15)} +
    labs(x="Predictor value",
         y="Relative average prediction") +
    theme_classic() +
    theme(axis.text=element_text(size=6),
          strip.text=element_text(size=7))
  ggsave(glue("figs/pub/ALE_{y}.png"), p, width=14, height=14*1.41, dpi=500)
  
  p <- adp_df |>
    ggplot(aes(x, value, colour=varTypeClean)) + 
    geom_hline(yintercept=0, linewidth=0.4, colour="black") +
    geom_line(aes(group=name), alpha=0.5, linewidth=0.2) + 
    geom_line(data=adp_ens, linewidth=0.5, colour="black") +
    scale_colour_brewer(type="qual", palette=2, guide="none") +
    {if(y_i$type[i]=="hab") facet_wrap(~variable_ordered, scales="free_x", nrow=14)} +
    {if(y_i$type[i]=="tox") facet_wrap(~variable_ordered, scales="free_x", nrow=15)} +
    scale_y_continuous(limits=range(adp_ens$value), oob=scales::oob_keep) +
    labs(x="Predictor value",
         y="Relative average prediction") +
    theme_classic() +
    theme(axis.text=element_text(size=6),
          strip.text=element_text(size=7))
  ggsave(glue("figs/pub/ALE_{y}_ensLim.png"), p, width=14, height=14*1.41, dpi=500)
  
  })
}

