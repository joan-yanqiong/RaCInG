#' Compute weights for ligand-receptor pairs
#'
#' @description Patient-specific LR pair weights were defined as the minimum of
#' log2(TPM+1) expression 613 of the ligand and the receptor, hypothesizing
#' that the expression of the gene at the lower level 614 limits the LR
#' binding affinity.
#' @param i Index of ligand-receptor pair
#' @param log2_tpm Log2-transformed TPM matrix
#' @param lr_pairs Data frame containing ligand-receptor pairs
#' @return Data frame containing weights for each patient
#' @export
#'
#' @examples  dontrun{
#' compute_lr_weight(1, log2_tpm, lr_pairs)}
compute_lr_weight <- function(i, log2_tpm, lr_pairs) {
    pos_lr <- match(c(lr_pairs[i, "Ligand"],
    lr_pairs[i, "Receptor"]), rownames(log2_tpm))
    if (sum(is.na(pos_lr)) > 0) {
        by_patient <- t(data.frame(rep(NA, ncol(log2_tpm))))
    } else {
        # When a ligand or receptor is not found, NA value should be returned.
        by_patient <- t(data.frame(apply(log2_tpm[pos_lr, ], 2, min)))
    }
    rownames(by_patient) <- lr_pairs[i, "interaction"]
    return(by_patient)
}