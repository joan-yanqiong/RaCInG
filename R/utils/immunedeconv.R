cellfrac_quantiseq %<-% immunedeconv::deconvolute_quantiseq(RNA_tpm,
    tumor = TRUE,
    arrays = FALSE,
    scale_mrna = TRUE
)
cellfrac_mcpcounter %<-% immunedeconv::deconvolute_mcp_counter(RNA_tpm,
    feature_types = "HUGO_symbols"
)

cellfrac_xcell %<-% immunedeconv::deconvolute_xcell(RNA_tpm, arrays = FALSE)

# cellfrac_epic %<-% immunedeconv::deconvolute_epic(RNA_tpm,
#     tumor = TRUE,
#     scale_mrna = TRUE
# )

