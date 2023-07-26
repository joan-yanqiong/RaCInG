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
    args$input_dir <- glue("{here::here()}/output/20_deconv")
    args$output_dir <- glue("{here::here()}/output/21_deconv_post")
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

# ---- Setup ----
cellfrac_quantiseq <- readRDS(glue("{args$input_dir}/{args$cancer_type}_cellfrac_quantiseq.rds"))
cellfrac_mcpcounter <- readRDS(glue("{args$input_dir}/{args$cancer_type}_cellfrac_mcpcounter.rds"))
cellfrac_xcell <- readRDS(glue("{args$input_dir}/{args$cancer_type}_cellfrac_xcell.rds"))
cellfrac_epic <- readRDS(glue("{args$input_dir}/{args$cancer_type}_cellfrac_epic.rds"))
cellfrac_timer <- readRDS(glue("{args$input_dir}/{args$cancer_type}_cellfrac_timer.rds"))

# ---- Processing  ----
# Determining consensus score for DCs
# 1. Scale to [0, 1]
# 2. Compute median
# 3. Scale to range of quantiseq

# Correct for quantiseq DCs using xCell
quantiseq_dc <- cellfrac_quantiseq["DCs", ]
quantiseq_dc[which(cellfrac_xcell["DCs", ] == 0)] <- 0

# 2. Compute median
dc_median <-
    apply(cbind(
        range01(cellfrac_timer),
        range01(cellfrac_mcpcounter),
        range01(cellfrac_xcell["DCs", ])
    ), 1, median)

# 3. Scale to range of quantiseq
dc_range <- range(quantiseq_dc)
dc <- dc_median * (dc_range[2] - dc_range[1]) + dc_range[1]

cell_fractions <- data.frame(rbind(cellfrac_quantiseq[rownames(cellfrac_quantiseq) != "DCs", ], cellfrac_epic, dc))
rownames(cell_fractions)[rownames(cell_fractions) == "dc"] <- "DCs"

# Combine M1 and M2 macrophages = macrophages
cell_fractions["Macrophages", ] <- cell_fractions["M1 macrophages", ] + cell_fractions["M2 macrophages", ]
cell_fractions <- cell_fractions[!rownames(cell_fractions) %in% c("M1 macrophages", "M2 macrophages"), ]

# Correct fractions using xCell
cell_fractions["Tregs", ] <- ifelse(cellfrac_xcell["Tregs", ] == 0, 0, cell_fractions["Tregs", ])
cell_fractions["NK cells", ] <- ifelse(cellfrac_xcell["NK cells", ] == 0, 0,
    cell_fractions["NK cells", ]
)

# Re-scaling the data to sum to 1
cell_fractions_norm <- cell_fractions
cell_fractions_norm <- cell_fractions_norm / apply(cell_fractions, 1, sum)

cell_fractions_norm <- t(cell_fractions_norm)

colnames(cell_fractions_norm) <- c("CD8+ T", "B", "Treg", "CAF", "NK", "Tumor", "Endo", "DC", "M")

# ---- Save results ----
log_info("Saving results...")
saveRDS(cell_fractions_norm, glue("{args$output_dir}/{args$cancer_type}_TMEmod_cell_fractions.rds"))

write.csv(cell_fractions_norm, glue("{args$output_dir}/{args$cancer_type}_TMEmod_cell_fractions.csv"))
