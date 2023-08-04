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
    parser$add_argument("-i", "--input_dir",
        type = "character",
        default = NULL, help = "Path to directory"
    )
    parser$add_argument("-c", "--cancer_type",
        type = "character",
        default = NULL, help = "Cancer type"
    )
    args <- parser$parse_args()
} else {
    # Provide arguments here for local runs
    args <- list()
    args$log_level <- 5
    args$input_dir <- glue("{here::here()}/output/GBM/40_immune_response")
    args$output_dir <- glue("{here::here()}/output/GBM/40_various_scores")
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
pacman::p_load(easier)

RNA_tpm <- readRDS(glue("{here::here()}/output/GBM/10_preprocessing/GBM_tpm.rds"))
RNA_counts <- readRDS(glue("{here::here()}/output/GBM/10_preprocessing/GBM_counts.rds"))


# Compute immune response scores
hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response_scores <- compute_scores_immune_response(
    RNA_tpm = RNA_tpm,
    # selected_scores = hallmarks_of_immune_response
)
saveRDS(immune_response_scores, glue("{args$output_dir}/{args$cancer_type}_immune_response_scores.rds"))


cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
saveRDS(cell_fractions, glue("{args$output_dir}/{args$cancer_type}_cell_fractions.rds"))

# Compute pathway activities
pathway_activities <- compute_pathway_activity(
    RNA_counts = RNA_counts,
    remove_sig_genes_immune_response = TRUE
)
saveRDS(pathway_activities, glue("{args$output_dir}/{args$cancer_type}_pathway_activities.rds"))


# Compute TF activities
tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)

saveRDS(tf_activities, glue("{args$output_dir}/{args$cancer_type}_tf_activities.rds"))

# Cell-cell interaction scores
lrpair_weights <- compute_LR_pairs(
    RNA_tpm = RNA_tpm,
    cancer_type = "pancan"
)
saveRDS(lrpair_weights, glue("{args$output_dir}/{args$cancer_type}_lrpair_weights.rds"))
ccpair_scores <- compute_CC_pairs(
    lrpairs = lrpair_weights,
    cancer_type = "pancan"
)
saveRDS(ccpair_scores, glue("{args$output_dir}/{args$cancer_type}_ccpair_scores.rds"))

predictions <- predict_immune_response(
    pathways = pathway_activities,
    immunecells = cell_fractions,
    tfs = tf_activities,
    lrpairs = lrpair_weights,
    ccpairs = ccpair_scores,
    cancer_type = args$cancer_type,
    verbose = TRUE
)
saveRDS(predictions, glue("{args$output_dir}/{args$cancer_type}_predictions.rds"))
