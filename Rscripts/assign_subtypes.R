pacman::p_load(Seurat, tidyverse, readxl, glue, ggplot2)

tpm_path <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/RaCInG/output/GBM/10_preprocessing/GBM_tpm.rds"
counts_path <- "/Users/joankant/Desktop/gaitigroup/Users/Joan/RaCInG/output/GBM/10_preprocessing/GBM_counts.rds"
tpm <- readRDS(tpm_path)
counts <- readRDS(counts_path)


suppfile <- glue("{here::here()}/Data/Neftel2019_subtypes.xlsx")
signatures <- read_excel(path = suppfile, skip = 4)
signatures <- signatures %>% select(-c("G1/S", "G2/M"))

tpm <- log2(1 + tpm / 10)
# score_matrix <- scalop::sigScores(
#     m = log_expression_matrix,
#     sigs = MES_functional_list
# )

seurat_obj <- CreateSeuratObject(counts)
seurat_obj <- SetAssayData(object = seurat_obj, assay = "RNA", slot = "data", new.data = as.matrix(tpm))


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6703186/
list_of_signatures <- lapply(colnames(signatures), function(signature) {
    genes_with_garbage <- signatures %>% pull(signature)
    return(genes_with_garbage[!is.na(genes_with_garbage)])
})

names(list_of_signatures) <- colnames(signatures)

seurat_obj <- AddModuleScore(object = seurat_obj, features = list_of_signatures, name = "Neftel2019_subtypes", assay = "RNA")

signature_scores <- seurat_obj@meta.data %>% select(contains("Neftel"))
colnames(signature_scores) <- colnames(signatures)

signature_scores <- signature_scores %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(
        cols = colnames(signature_scores),
        names_to = "signature", values_to = "score"
    ) %>%
    # Get the subtype with the highest score for each sample
    group_by(sample_id) %>%
    top_n(1, score) %>%
    rowwise() %>%
    # Adapt labelling: 3 categories, 4 categories
    mutate(
        Neftel3 = case_when(
            signature == "MES1" | signature == "MES2" ~ "MES-like",
            signature == "NPC1" | signature == "NPC2" | signature == "OPC" ~ "NPC & OPC-like",
            signature == "AC" ~ "AC-like",
        ),
        Neftel4 = case_when(
            signature == "MES1" | signature == "MES2" ~ "MES-like",
            signature == "NPC1" | signature == "NPC2" ~ "NPC-like",
            signature == "AC" ~ "AC-like",
            signature == "OPC" ~ "OPC-like"
        )
    )

subtyping_long <- signature_scores %>% pivot_longer(cols = c("signature", "Neftel3", "Neftel4"), names_to = "method", values_to = "subtype")
g <- ggplot(data = subtyping_long) +
    facet_wrap(~method, scales = "free") +
    geom_bar(aes(x = subtype))

# signature_scores$Patient <- substr(signature_scores$sample_id, 1, 12)

signature_scores <- readRDS("/Users/joankant/Desktop/gaitigroup/Users/Joan/RaCInG/output/GBM/GBM_Neftel2019_subtypes.rds")

test <- signature_scores %>%
    mutate(Patient = str_replace_all(substr(sample_id, 1, 15), "\\.", "-"), sample_id = str_replace_all(sample_id, "\\.", "-"))
test <- test[, c("Patient", "Neftel3", "Neftel4", "signature", "score", "sample_id")]

ggsave(glue("{here::here()}/output/GBM/figures/TCGA_GBM_Neftel2019_subtypes.pdf"), plot = g)

saveRDS(signature_scores, glue("{here::here()}/output/GBM/GBM_Neftel2019_subtypes.rds"))

write.csv(test, glue("{here::here()}/output/GBM/GBM_Neftel2019_subtypes.csv"), row.names = FALSE)
