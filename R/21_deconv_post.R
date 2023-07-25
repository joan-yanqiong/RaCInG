if (!require("pacman")) install.packages("pacman")

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(argparse, glue, data.table, tidyverse, stringr)

# Parse user input
parser <- ArgumentParser(description = "Deconvolution Post-Processing")
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

# Load additional libraries (use pacman::p_load() or pacman::p_load_gh() to load
# libraries)

create_dir(args$output_dir)

# ---- Setup ----
cellfrac_quantiseq <- readRDS(glue("{args$output_dir}/{args$cancer_type}_cellfrac_quantiseq.rds"))
cellfrac_mcpcounter <- readRDS(glue("{args$output_dir}/{args$cancer_type}_cellfrac_mcpcounter.rds"))
cellfrac_xcell <- readRDS(glue("{args$output_dir}/{args$cancer_type}_cellfrac_xcell.rds"))
cellfrac_epic <- readRDS(glue("{args$output_dir}/{args$cancer_type}_cellfrac_epic.rds"))
cellfrac_timer <- readRDS(glue("{args$output_dir}/{args$cancer_type}_cellfrac_timer.rds"))

# \' Scale values to [0, 1]
# \'
# \' @param x numeric vector
# \' @return vector with values scaled to [0, 1]
# \' @export
# \'
# \' @examples range01(c(1, 2, 3))

range01 <- function(x) {
    return(x - min(x)) / (max(x) - min(x))
}

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
