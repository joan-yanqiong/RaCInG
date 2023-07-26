#' Detect duplicate genes in expression matrix
#'
#' @param genes vector of gene symbols
#' @return vector of duplicate genes
#' @export
#'
#' @examples detect_duplicate_genes(c("A", "B", "A", "C"))
detect_duplicate_genes <- function(genes) {
    n_occurences <- table(genes)
    dupl_genes <- names(n_occurences[n_occurences > 1])
    return(dupl_genes)
}

#' Average expression of duplicate genes
#'
#' @param gene_expr gene expression matrix
#' @param dupl_genes vector of duplicate genes
#' @return gene expression matrix with averaged expression of duplicate genes
#' @export
#'
#' @examples avg_expr(gene_expr, c("A", "B", "C")))
#' @importFrom dplyr group_by summarise_all filter
avg_expr <- function(gene_expr) {
    return(gene_expr %>%
        filter(symbol %in% dupl_genes) %>%
        group_by(symbol) %>%
        summarise_all(mean))
}

#' Remove duplicate genes from expression matrix by averaging their expression
#'
#' @param gene_expr gene expression matrix
#' @return gene expression matrix without duplicate genes
#' @export
#'
#' @examples remove_dupl_genes(gene_expr)
#' @importFrom dplyr filter bind_rows
#' @importFrom tibble column_to_rownames
remove_dupl_genes <- function(gene_expr) {
    dupl_genes <- detect_duplicate_genes(gene_expr$symbol)
    if (length(dupl_genes) == 0) {
        return(gene_expr)
    }
    avg_expr_dupl <- avg_expr(gene_expr %>% filter(symbol %in% dupl_genes))
    return(gene_expr %>%
        filter(!symbol %in% dupl_genes) %>% bind_rows(avg_expr_dupl) %>% column_to_rownames("symbol"))
}
