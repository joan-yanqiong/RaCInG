pacman::p_load(tidyverse)
path <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/RaCInG/output/GBM/kernel/GBM_D.out"
lines <- readLines(path)
lines <- noquote(lines)
lines <- stringr::str_replace_all(lines, c("\\[" = "", "\\]" = "", "'" = ""))

meta <- lines[1]
patient_id <- NULL
celltypes <- NULL
extract_comm_proba <- NULL

track_celltypes <- list()
comm_scores <- list()
for (line in lines[2:length(lines)]) {
    extract_cell_types <- strsplit(line, " ")[[1]]
    extract_comm_proba <- strsplit(line, ",")[[1]]

    if (length(extract_cell_types) > 1) {
        celltypes <- extract_cell_types
    } else if (length(extract_comm_proba) > 1) {
        comm_proba <- extract_comm_proba
    } else {
        patient_id <- as.integer(line)
    }
    track_celltypes <- append(track_celltypes, data.frame(c(patient_id, celltypes)))
    comm_scores <- append(comm_scores, data.frame(as.numeric(c(patient_id, comm_proba))))
}

celltypes <- data.frame(do.call(rbind, track_celltypes), row.names = NULL) %>% distinct()

celltypes_oi <- c("B", "CAF", "CD8+ T", "DC", "Endo", "M", "NK", "Treg", "Tumor")
names(celltypes_oi) <- seq_len(length(celltypes_oi)) - 1

comm_scores_df <- data.frame(do.call(rbind, comm_scores), row.names = NULL)
colnames(comm_scores_df) <- c("patient_id", "source", "target", "comm_proba")
comm_scores_df$source <- unlist(sapply(comm_scores_df$source, function(id) {
    return((celltypes_oi[[as.character(id)]]))
}))
comm_scores_df$target <- unlist(sapply(comm_scores_df$target, function(id) {
    return((celltypes_oi[[as.character(id)]]))
}))
comm_scores_df <- data.frame(comm_scores_df %>% filter(!is.na(patient_id)))


# add metadata
subtypes <- readRDS("/Users/joankant/Desktop/gaitigroup/Users/Joan/RaCInG/output/GBM/GBM_Neftel2019_subtypes.rds")

subtypes_vec <- subtypes$Neftel3
names(subtypes_vec) <- as.character(seq_len(length(subtypes_vec)) - 1)


# comm_scores_df <- str_replace_all(comm_scores_df, subtypes_vec)
# comm_scores_df$subtype <-


comm_scores_df$subtype <- unlist(sapply(comm_scores_df$patient_id, function(id) {
    return((subtypes_vec[[as.character(id)]]))
}))

avg_comm <- comm_scores_df %>%
    group_by(subtype, source, target) %>%
    summarise(mean_comm_proba = mean(comm_proba))


avg_comm %>%
    group_by(subtype) %>%
    chordDiagram()

par(mfrow = c(1, 3))

# chordDiagram(avg_comm %>% filter(subtype == "AC-like") %>% ungroup() %>% select(-subtype))
# circos.clear()
# chordDiagram(avg_comm %>% filter(subtype == "MES-like") %>% ungroup() %>% select(-subtype))
# circos.clear()
# chordDiagram(avg_comm %>% filter(subtype == "NPC & OPC-like") %>% ungroup() %>% select(-subtype))

ac_like <- avg_comm %>%
    filter(subtype == "AC-like") %>%
    ungroup() %>%
    select(-subtype) %>%
    data.frame()
colnames(ac_like) <- c("from", "to", "value")
chordDiagram(ac_like, gap.degree = 0.01, order = celltypes_oi)
