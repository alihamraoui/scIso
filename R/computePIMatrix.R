#' Compute PI (Percent Isoform) Matrix from Isoform Counts
#'
#' This function calculates per-cell isoform usage (PI) values from an isoform-level
#' count matrix stored in a Seurat object. For each gene, isoform counts are divided
#' by the total counts of all isoforms belonging to that gene within each cell.
#'
#' @param seurat_obj A Seurat object containing isoform-level counts in one assay.
#' @param assay Character string specifying the assay that contains isoform counts.
#'   Default is `"ISO"`.
#' @param slot Character string indicating which slot of the assay to use
#'   (e.g., `"counts"` or `"data"`). Default is `"counts"`.
#'
#' @details
#' The function expects rownames in the form `GeneName..TranscriptID`. The part
#' before the first double dot (`..`) is interpreted as the gene identifier.
#' PI values are computed per isoform and per cell as:
#' \deqn{PI = count(isoform) / sum(counts of all isoforms from the same gene)}
#'
#' Cells with zero total counts for a given gene return `NA` PI values for its isoforms.
#'
#' @return A matrix of PI values (isoforms - cells) with \code{NA} where
#'   total gene counts are zero.
#'
#' @export
#'
#' @examples
#' library(Seurat)
#' iso_names <- c("G1..Tx1", "G1..Tx2", "G2..Tx1", "G2..Tx2")
#' counts <- matrix(c(10, 5, 3, 7, 0, 4, 2, 6),
#'                  nrow = 4,
#'                  dimnames = list(iso_names, paste0("Cell", 1:2)))
#' seurat_obj <- CreateSeuratObject(counts = counts, assay = "ISO")
#'
#' pi <- computePIMatrix(seurat_obj)
#' head(pi)
#'
computePIMatrix <- function(seurat_obj,
                              assay = "ISO",
                              slot = "counts") {

  # --- Basic checks ---
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }

  iso_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  if (nrow(iso_mat) == 0L) {
    stop("Isoform matrix is empty.")
  }

  base::message("[computePIMatrix] Extracted assay '", assay, "' (",
                base::nrow(iso_mat), " isoforms - ", base::ncol(iso_mat), " cells).")

  # --- Compute PSI ---
  gene_names <- base::sub("\\.\\..*$", "", base::rownames(iso_mat))
  gene_split <- base::split(base::seq_along(gene_names), gene_names)
  gene_totals <- base::do.call(base::rbind, base::lapply(gene_split, function(idx) {
    Matrix::colSums(iso_mat[idx, , drop = FALSE])
  }))
  gene_totals_mat <- gene_totals[gene_names, , drop = FALSE]
  pi_mat <- iso_mat / gene_totals_mat

  base::message("[computePIMatrix] PI matrix computed successfully (dim: ",
                base::nrow(pi_mat), " - ", base::ncol(pi_mat), ").")

  #pi_mat[is.na(pi_mat)] <- 0
  return(as.matrix(pi_mat))
}
