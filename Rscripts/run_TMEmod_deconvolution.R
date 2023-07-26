# # --------------------------------------------------- #
# #  Script to run deconvolution computational methods   #
# #  Francesca Finotello wrote this code.                #
# # --------------------------------------------------- #
# pacman::p_load(singscore)
# pacman::p_load_gh("cansysbio/ConsensusTME@v0.0.1.9000", "omnideconv/immunedeconv@2.1.0", "dviraran/xCell@1.1.0")

# run_TMEmod_deconvolution <- function(RNA_tpm, cancer_type) {
#   cellfrac_consensusTME <- ConsensusTME::consensusTMEAnalysis(as.matrix(RNA_tpm),
#     cancer = cancer_type,
#     statMethod = "ssgsea"
#   )

#   cellfrac_quanTIseq <- immunedeconv::deconvolute_quantiseq(RNA_tpm,
#     tumor = TRUE,
#     arrays = FALSE,
#     scale_mrna = TRUE
#   )

#   cellfrac_EPIC <- immunedeconv::deconvolute_epic(RNA_tpm,
#     tumor = TRUE,
#     scale_mrna = TRUE
#   )

#   cellfrac_mcpcounter <- immunedeconv::deconvolute_mcp_counter(RNA_tpm,
#     feature_types = "HUGO_symbols"
#   )

#   cellfrac_xCell <- immunedeconv::deconvolute_xcell(RNA_tpm,
#     arrays = FALSE
#   )

#   cellfrac_TIMER <- immunedeconv::deconvolute_timer(RNA_tpm,
#     indications = rep(cancer_type, ncol(RNA_tpm))
#   )

#   return(list(
#     quanTIseq = cellfrac_quanTIseq,
#     EPIC = cellfrac_EPIC,
#     mcpcounter = cellfrac_mcpcounter,
#     xCell = cellfrac_xCell,
#     TIMER = cellfrac_TIMER,
#     consensusTME = cellfrac_consensusTME
#   ))
# }
