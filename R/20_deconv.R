if (!require("pacman")) install.packages("pacman")

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(argparse, glue, data.table, tidyverse, stringr)

# Parse user input
parser <- ArgumentParser(description = "Deconvolution of cell types")
parser$add_argument("-ll", "--log_level",
    type = "integer",
    default = "4", help = "Log level: 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, 5=DEBUG"
)
parser$add_argument("-o", "--output_dir",
    type = "character",
    default = NULL, help = "Directory to save output"
)

args <- parser$parse_args()

if (interactive()) {
    # 	Provide arguments here for local runs
    print("Running interactively...")
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/20_deconv")
    args$input_file <- glue("{here::here()}/output/BLCA_tpm.rds")
    args$cancer_type <- "BLCA"
} else {
    # 	Set temporary working directory based on path of script
    print("Running from command line/terminal...")
    setwd(dirname(str_split(commandArgs(trailingOnly = FALSE)
    [grep("--file=", commandArgs(trailingOnly = FALSE))], "=", simplify = TRUE)[2]))
}
pacman::p_load(here)
# Set final working directory
setwd(here::here())
# Load standard user-defined functions
source(glue("{here::here()}/R/utils/utils.R"))
source(glue("{here::here()}/R/utils/run_TMEmod_deconvolution.R"))
# Set up logging
logr <- init_logging(log_level = args$log_level)

# Load additional libraries (use pacman::p_load() or pacman::p_load_gh() to load
# libraries)
pacman::p_load(parallel)
pacman::p_load_gh("omnideconv/immunedeconv")
# plan(multisession)

# ---- Constants ----
quantiseq_oi <- c("T.cells.CD8", "B.cells", "Tregs", "Macrophages.M1", "Macrophages.M2", "Dendritic.cells")
epic_oi <- c("CAFs", "NKcells", "otherCells", "Endothelial")
# NKcells from epic_oi and Tregs from quantiseq set to zero if xCell is zero
xcell_oi <- c("NK cells", "Tregs", "DC")
mcp_counter <- c("Myeloid dendritic cells")
timer_oi <- c("DC")


# Read RNA-seq data
RNA_tpm_df <- readRDS(args$input_file)
RNA_tpm <- data.matrix(RNA_tpm_df %>% select(-symbol))
rownames(RNA_tpm) <- RNA_tpm_df$symbol

# ---- Compute cell abundances ----
log_info("Compute cell abundances using")


cellfrac_quantiseq <- immunedeconv::deconvolute_quantiseq(RNA_tpm,
    tumor = TRUE,
    arrays = FALSE,
    scale_mrna = TRUE
)
cellfrac_mcpcounter <- immunedeconv::deconvolute_mcp_counter(RNA_tpm,
    feature_types = "HUGO_symbols"
)

cellfrac_xcell <- immunedeconv::deconvolute_xcell(RNA_tpm, arrays = FALSE)

cellfrac_epic <- immunedeconv::deconvolute_epic(RNA_tpm,
    tumor = TRUE,
    scale_mrna = TRUE
)

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
