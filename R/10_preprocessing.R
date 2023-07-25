if (!require("pacman")) install.packages("pacman")

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(argparse, glue, data.table, tidyverse, stringr)

# Parse user input
parser <- ArgumentParser(description = "Pre-processing")
parser$add_argument("-ll", "--log_level",
    type = "integer",
    default = "4", help = "Log level: 1=FATAL, 2=ERROR, 3=WARN, 4=INFO, 5=DEBUG"
)
parser$add_argument("-o", "--output_dir",
    type = "character",
    default = NULL, help = "Directory to save output"
)

parser$add_argument("-i", "--input_file",
    type = "character",
    default = NULL, help = "Input file (txt)"
)

args <- parser$parse_args()

if (interactive()) {
    # 	Provide arguments here for local runs
    print("Running interactively...")
    args <- list()
    args$log_level <- 5
    args$output_dir <- glue("{here::here()}/output/10_preprocessing")
    args$input_file <- glue("{here::here()}/Data/BLCA_rnaseq.txt")
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
pacman::p_load_gh("maialab/hgnc")
source(glue("{here::here()}/R/utils/preprocessing.R"))

# ---- Initialisation ----
create_dir(args$output_dir)
file <- fread(args$input_file)
hgnc_dataset <- import_hgnc_dataset()

# ---- Pre-processing ----
# Formatting
genes <- rownames(as.matrix(file, rownames = 1))
symbol <- str_split(genes, "\\|", simplify = TRUE)[, 1]
hgnc_id2 <- str_split(genes, "\\|", simplify = TRUE)[, 2]

colnames(file) <- paste(colnames(file), file[1, ], sep = "__")
file <- data.frame(file)
file$symbol <- symbol
# file <- file[-c(1), ]

# Convert to long format
file_long <- file %>%
    select(-c(Hybridization.REF__gene_id)) %>%
    reshape2::melt(id.vars = c("symbol"))

# For speeding up, not using separate()
tmp <- data.frame(str_split(file_long$variable, "__", simplify = TRUE))

# Combine and check for valid HGNC ids
file_long <- cbind(file_long, tmp) %>%
    filter(symbol %in% hgnc_dataset$symbol) %>%
    rename(sample_id = "X1", type = "X2", expr = "value") %>%
    select(-c(variable))

log_info("Compute TPM and convert to wide format...")
tpm <- file_long %>%
    filter(type == "scaled_estimate") %>%
    select(-type) %>%
    mutate(expr = as.numeric(expr) * 1e6) %>%
    ungroup() %>%
    reshape2::dcast(symbol ~ sample_id, value.var = "expr")

log_info("Convert to wide format...")
counts <- file_long %>%
    filter(type == "raw_count") %>%
    select(-type) %>%
    mutate(value = as.numeric(expr)) %>%
    ungroup() %>%
    reshape2::dcast(symbol ~ sample_id, value.var = "expr")

log_info("Compute mean of duplicate genes...")

tpm <- remove_dupl_genes(tpm) %>% column_to_rownames("symbol")
counts <- remove_dupl_genes(counts) %>% column_to_rownames("symbol")

log_info("Save gene expression data...")
saveRDS(counts, glue("{args$output_dir}/{get_name(args$cancer_type)}_counts.rds"))
saveRDS(tpm, glue("{args$output_dir}/{get_name(args$cancer_type)}_tpm.rds"))

log_info("COMPLETED!")
