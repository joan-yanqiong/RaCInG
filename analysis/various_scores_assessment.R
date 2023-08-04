pacman::p_load(here, glue, tidyverse, ggplot2, ggpubr)
pacman::p_load_gh("jokergoo/ComplexHeatmap")
subtypes <- readRDS(glue("{here::here()}/output/GBM/GBM_Neftel2019_subtypes.rds")) %>% as.data.frame()

output_dir <- glue("{here::here()}/output/GBM/40_various_scores/")

tf_activity <- readRDS(glue("{output_dir}/GBM_tf_activities.rds")) %>% as.data.frame()

pathway_activities <- readRDS(glue("{output_dir}/GBM_pathway_activities.rds")) %>% as.data.frame()

Heatmap(pathway_activities)
