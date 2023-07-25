if (!require("pacman")) install.packages("pacman")

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(argparse, glue, data.table, tidyverse, stringr)

# Parse user input
parser <- ArgumentParser(description = "Compute LR-pair weights")
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
    args$output_dir <- glue("{here::here()}/output/30_compute_lr_pair_weights")
    args$input_file <- glue("{here::here()}/output/BLCA_tpm.rds")
    args$ccc_table_path <- "/Users/joankant/Library/CloudStorage/OneDrive-UHN/RaCInG/Data/expressed_ccc_table.csv"
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
# Set up logging
logr <- init_logging(log_level = args$log_level)

# Load additional libraries (use pacman::p_load() or pacman::p_load_gh() to load libraries)

tpm <- readRDS(args$input_file)
log2_tpm <- tpm %>% mutate_all(function(x) log2(as.numeric(x) + 1))

# curated cell-cell communication table
ccc_table_filtered <- read.csv(file = args$ccc_table_path)

# Extract unique ligand-receptor pairs
lr_pairs <- ccc_table_filtered %>%
    mutate(interaction = paste0(Ligand, "_", Receptor)) %>%
    select(interaction, Ligand, Receptor) %>%
    distinct(interaction, .keep_all = TRUE)

# \' Compute weights for ligand-receptor pairs
# \' Patient-specific LR pair weights were defined as the minimum of the
# \ 'log2(TPM+1) expression 613 of the ligand and the receptor, hypothesizing
# \' that the expression of the gene at the lower level 614 limits the LR
# \' binding affinity.
# \' @param i Index of ligand-receptor pair
# \' @param log2_tpm Log2-transformed TPM matrix
# \' @param lr_pairs Data frame containing ligand-receptor pairs
# \' @return Data frame containing weights for each patient
# \' @export
# \'
# \' @examples   \dontrun{
# \' compute_lr_weight(1, log2_tpm, lr_pairs)
compute_lr_weight <- function(i, log2_tpm, lr_pairs) {
    pos_lr <- match(c(lr_pairs[i, "Ligand"], lr_pairs[i, "Receptor"]), rownames(log2_tpm))
    if (sum(is.na(pos_lr)) > 0) {
        by_patient <- t(data.frame(rep(NA, ncol(log2_tpm))))
    } else {
        # When a ligand or receptor is not found, NA value should be returned.
        by_patient <- t(data.frame(apply(log2_tpm[pos_lr, ], 2, min)))
    }
    rownames(by_patient) <- lr_pairs[i, "interaction"]
    return(by_patient)
}

lr_weights <- do.call(
    rbind,
    lapply(seq_len(nrow(lr_pairs)),
        compute_lr_weight,
        lr_pairs = lr_pairs, log2_tpm = log2_tpm
    )
)

colnames(lr_weights) <- str_replace_all(colnames(lr_weights), "\\.", "-")

lr_weights <- t(lr_weights)
saveRDS(lr_weights,
    file = glue("{args$output_dir}/{args$cancer_type}_LRpairs_weights_min.rds")
)
write.csv(lr_weights,
    file = glue("{args$output_dir}/{args$cancer_type}_LRpairs_weights_min.csv")
)

log_info("COMPLETED!")
