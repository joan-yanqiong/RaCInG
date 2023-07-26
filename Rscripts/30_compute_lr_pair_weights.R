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
        description = "Get metadata",
    )
    parser$add_argument("-i", "--input_file",
        type = "character",
        default = NULL, help = "Path to Seurat object"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/30_compute_lr_pair_weights")
    args$input_file <- glue("{here::here()}/output/BLCA_tpm.rds")
    args$ccc_table_path <- glue("{here::here()}/Data/expressed_ccc_table.csv")
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

tpm <- readRDS(args$input_file)
log2_tpm <- tpm %>% mutate_all(function(x) log2(as.numeric(x) + 1))

# curated cell-cell communication table
ccc_table_filtered <- read.csv(file = args$ccc_table_path)

# Extract unique ligand-receptor pairs
lr_pairs <- ccc_table_filtered %>%
    mutate(interaction = paste0(Ligand, "_", Receptor)) %>%
    select(interaction, Ligand, Receptor) %>%
    distinct(interaction, .keep_all = TRUE)

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
