# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Set working directory
cmd_args <- commandArgs(trailingOnly = FALSE)
has_script_filepath <- startsWith(cmd_args, "--file=")
if (sum(has_script_filepath)) {
    setwd(dirname(unlist(strsplit(cmd_args[has_script_filepath], "=")))[2])
}

# Load libraries
pacman::p_load(glue, data.table, tidyverse, stringr)
devtools::load_all("./", export_all = FALSE)

if (!interactive()) {
    # Define input arguments when running from bash
    parser <- setup_default_argparser(
        description = "Compute cell type abundances",
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to Seurat object"
    )
    parser$add_argument("-c", "--cancer_type", type = "character", default = NULL, help = "Cancer type")
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/20_deconv")
    args$input_file <- glue("{here::here()}/output/10_preprocessing/GBM_tpm.rds")
    args$cancer_type <- "GBM"
}

# Set up logging
logr <- init_logging(log_level = args$log_level)
log_info(ifelse(interactive(),
    "Running interactively...",
    "Running from command line/terminal..."
))

log_info("Create output directory...")
create_dir(args$output_dir)

# Load additional libraries
pacman::p_load(parallel)
pacman::p_load_gh("omnideconv/immunedeconv")

# ---- Constants ----
quantiseq_oi <- c("T.cells.CD8", "B.cells", "Tregs", "Macrophages.M1", "Macrophages.M2", "Dendritic.cells")
epic_oi <- c("CAFs", "NKcells", "otherCells", "Endothelial")
# NKcells from epic_oi and Tregs from quantiseq set to zero if xCell is zero
xcell_oi <- c("NK cells", "Tregs", "DC")
mcp_counter <- c("Myeloid dendritic cells")
timer_oi <- c("DC")


# Read RNA-seq data
RNA_tpm_df <- readRDS(args$input_file)
RNA_tpm <- data.matrix(RNA_tpm_df)

# ---- Compute cell abundances ----
log_info("Compute cell abundances using")

log_info("Quantiseq...")
cellfrac_quantiseq <- immunedeconv::deconvolute_quantiseq(RNA_tpm,
    tumor = TRUE,
    arrays = FALSE,
    scale_mrna = TRUE
)
log_info("MCP-counter...")
cellfrac_mcpcounter <- immunedeconv::deconvolute_mcp_counter(RNA_tpm,
    feature_types = "HUGO_symbols"
)
log_info("xCell...")
cellfrac_xcell <- immunedeconv::deconvolute_xcell(RNA_tpm, arrays = FALSE)
log_info("EPIC...")
cellfrac_epic <- immunedeconv::deconvolute_epic(RNA_tpm,
    tumor = TRUE,
    scale_mrna = TRUE
)
log_info("TIMER...")
cellfrac_timer <- immunedeconv::deconvolute_timer(RNA_tpm,
    indications = rep(args$cancer_type, ncol(RNA_tpm))
)


# ---- Select cell types for each method ----
cellfrac_quantiseq <- cellfrac_quantiseq[quantiseq_oi, ]
cellfrac_epic <- cellfrac_epic[epic_oi, ]
cellfrac_mcpcounter <- cellfrac_mcpcounter[mcp_counter, ]
cellfrac_xcell <- cellfrac_xcell[xcell_oi, ]
cellfrac_timer <- cellfrac_timer[timer_oi, ]

rownames(cellfrac_quantiseq) <- c("CD8 T cells", "B cells", "Tregs", "M1 macrophages", "M2 macrophages", "DCs")
rownames(cellfrac_epic) <- c("CAFs", "NK cells", "Tumor cells", "Endothelial cells")
rownames(cellfrac_xcell) <- c("NK cells", "Tregs", "DCs")

# TODO remove after finishing
saveRDS(cellfrac_quantiseq, glue("{args$output_dir}/{args$cancer_type}_cellfrac_quantiseq.rds"))
saveRDS(cellfrac_mcpcounter, glue("{args$output_dir}/{args$cancer_type}_cellfrac_mcpcounter.rds"))
saveRDS(cellfrac_xcell, glue("{args$output_dir}/{args$cancer_type}_cellfrac_xcell.rds"))
saveRDS(cellfrac_epic, glue("{args$output_dir}/{args$cancer_type}_cellfrac_epic.rds"))
saveRDS(cellfrac_timer, glue("{args$output_dir}/{args$cancer_type}_cellfrac_timer.rds"))
