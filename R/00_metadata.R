if (!require("pacman")) install.packages("pacman")

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(argparse, glue, data.table, tidyverse, stringr)

# Parse user input
parser <- ArgumentParser(description = "Metadata")
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
    args$output_dir <- glue("{here::here()}/output/")
    args$input_file <- "/Users/joankant/Library/CloudStorage/OneDrive-UHN/TCGA-GBM/data/clinical/GBM.clin.merged.picked.txt"
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

meta <- data.frame(fread(args$input_file, sep = "\t", header = FALSE))
rownames(meta) <- meta$V1
meta <- t(meta)
rownames(meta) <- NULL
meta <- meta[-c(1), ]
meta <- data.frame(meta)
hist_types <- unique(meta$histological_type)
