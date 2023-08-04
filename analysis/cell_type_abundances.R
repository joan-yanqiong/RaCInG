pacman::p_load(here, glue, tidyverse, ggplot2, ggpubr)
subtypes <- readRDS(glue("{here::here()}/output/GBM/GBM_Neftel2019_subtypes.rds")) %>% as.data.frame()

cellfractions <- readRDS(glue("{here::here()}/output/GBM/21_deconv_post/GBM_TMEmod_cell_fractions.rds")) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id")

cellfractions_long <- cellfractions %>% pivot_longer(values_to = "fraction", cols = -sample_id, names_to = "cell type")


celltype_data <- cellfractions_long %>% left_join(subtypes, by = "sample_id", )

expand.grid(c("AC-like", "NPC & OPC-like", "MES-like"), c("AC-like", "NPC & OPC-like", "MES-like"))

comparisons <- list(c("AC-like", "NPC & OPC-like"), c())

comps <- t(combn(c("AC-like", "NPC & OPC-like", "MES-like"), 2, simplify = TRUE))

comps <- lapply(seq_len(nrow(comps)), function(x) {
    return(c(comps[x, ]))
})

p <- ggboxplot(
    data = celltype_data, x = "Neftel3", y = "fraction",
    color = "Neftel3", facet.by = "cell type", add = "jitter",
    scales = "free"
) + stat_compare_means(comparisons = comps, method = "wilcox.test", label = "p.signif")
print(p)

ggsave(plot = p, filename = glue("{here::here()}/output/figures/GBM_TMEmod_cell_fractions_Neftel3.pdf"), width = 12, height = 14)

subtypes %>%
    select(sample_id, Neftel3) %>%
    group_by(Neftel3) %>%
    summarise(n = n())

nsamples <- ggplot(data = subtypes, aes(x = Neftel3)) +
    geom_bar(stat = "count") +
    stat_count(
        geom = "text", colour = "white", size = 10,
        aes(label = after_stat(count)), position = position_stack(vjust = .8)
    )
