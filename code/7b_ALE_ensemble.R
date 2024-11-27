# Ensemble ALE plots

library(tidyverse)
library(tidymodels)
library(ingredients)
library(sevcheck)
library(glue)
source("code/00_fn.R")

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


for(i in seq_along(y_i)) {
  y <- y_i$abbr[i]
  
  if(file.exists(glue("out/0_init/compiled/{y}_adp_Ensemble.rds"))) next
  cat("Starting", y, "\n")
  # For each monitoring target, load and merge all CP profiles
  # For variables not appearing in a candidate model, use constant (mean pr).
  cp_df <- dirrf(glue("{out.dir}/compiled/"), glue("{y}_cp.[^data]")) |>
    map(readRDS) |>
    reduce(full_join, by = join_by(`_vname_`, `_ids_`, cp_x_id)) |>
    mutate(obsid=row_number()) |>
    mutate(across(ends_with("A1"), ~replace_na(.x, mean(.x, na.rm=T))))
  cat("Loaded CP\n")
  cp_attr <- readRDS(glue("{out.dir}/compiled/16-Avg1_Xf1_XN1_Del1/{y}_cp_data_attr.rds")) |>
    select(y, siteid, date, `_ids_`) |>
    mutate(`_ids_`=as.character(`_ids_`))
  
  
  # Calculate ensemble from candidate CP probabilities
  cat("Calculating ensemble\n")
  cp_df_new <- calc_ensemble(out.ls=list(alert=cp_df |> inner_join(cp_attr)),
                             wt.ls=NULL,
                             resp="alert",
                             y_i.i=y_i[i,],
                             method="GLM_oos",
                             out.path=glue("{out.dir}/ensembles/"))
  
  # Merge with CP data
  cp_df <- inner_join(readRDS(glue("{out.dir}/compiled/16-Avg1_Xf1_XN1_Del1/{y}_cp_data.rds")), 
                      cp_df_new |> select(`_vname_`, `_ids_`, cp_x_id, ensGLM2_alert_A1)) |>
    rename(`_yhat_`=ensGLM2_alert_A1) |>
    mutate(`_label_`="Ensemble")
  
  attr(cp_df, "observations") <- readRDS(glue("{out.dir}/compiled/1-Avg0_Xf0_XN0_Del0/{y}_cp_data_attr.rds"))
  class(cp_df) <-  c("ceteris_paribus_explainer", "tbl_df", "tbl", "data.frame")
  adp <- aggregate_profiles(cp_df, variables = unique(cp_df$`_vname_`), type = "accumulated", variable_type = "numerical")
  adp_ens <- adp |>
    set_names(c("variable", "model", "x", "yhat", "id")) |>
    select(variable, x, yhat)
  saveRDS(adp_ens, glue("out/0_init/compiled/{y}_adp_Ensemble.rds"))
}

