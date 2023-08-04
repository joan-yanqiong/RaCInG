pacman::p_load(here, glue, tidyverse, ggplot2, ggpubr)
subtypes <- readRDS(glue("{here::here()}/output/GBM/GBM_Neftel2019_subtypes.rds")) %>% as.data.frame()

subtypes %>%
    select(sample_id, Neftel3) %>%
    group_by(Neftel3) %>%
    summarise(n = n())

nsamples <- ggplot(data = subtypes, aes(x = Neftel3, fill = Neftel3)) +
    geom_bar(stat = "count") +
    stat_count(
        geom = "text", colour = "white", size = 10,
        aes(label = after_stat(count)), position = position_stack(vjust = .8)
    )

ggsave(plot = nsamples, filename = glue("{here::here()}/output/figures/GBM_Neftel2019_subtypes.pdf"), width = 12, height = 14)
