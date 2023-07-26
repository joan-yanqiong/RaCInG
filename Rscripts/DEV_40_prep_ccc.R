if (!require("pacman")) install.packages("pacman")

# Unload all previously loaded packages + remove previous environment
rm(list = ls(all = TRUE))
pacman::p_unload()

# Load libraries
pacman::p_load(argparse, glue, data.table, tidyverse, stringr)

# Parse user input
parser <- ArgumentParser(description = "Prepare for CCC analysis")
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
pacman::p_load(reshape2)
pacman::p_load_gh("saezlab/OmnipathR", "saezlab/liana", "olapuentesantana/easier")

# Obtain curated CCC interactions
# Retrieve consensus resource, curated  (in the context of CCC) ligand-receptor interactions
curated_omni <- liana::get_curated_omni()
# Remove interactions from ECM
curated_omni <- curated_omni %>%
    filter(category_intercell_source %in% c("cell_surface_ligand", "ligand") &
        category_intercell_target %in% c("receptor", "adhesion")) %>%
    ## Decomplexify (or split all complexes into subunits)
    liana::decomplexify(columns = c("source_genesymbol", "target_genesymbol"))

#  Ligand-receptor pairs
lr_pairs_curated <- curated_omni %>%
    mutate(interaction_pair = paste0(source_genesymbol, "_", target_genesymbol)) %>%
    select(interaction_pair) %>%
    unique() %>%
    pull()

# Final knowledge data.frame
intercell_knowledge <- curated_omni
