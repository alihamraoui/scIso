#' Boxplot of isoform expression across conditions
#'
#' This function visualizes isoform expression (or PI if provided upstream)
#' across a categorical variable (e.g. cell cycle phase), optionally faceted by
#' another metadata column, and performs per-panel statistical comparison
#' between isoforms using Wilcoxon tests.
#'
#' @param seurat_obj A Seurat object.
#' @param isoforms Character vector of isoform feature names to include
#'   (e.g. c("Cdkn2c..ENSMUST00000097921", "Cdkn2c..ENSMUST00000063531")).
#' @param group_col Metadata column used on the x-axis (e.g. "Phase").
#' @param facet_col Optional metadata column used to facet the plot
#'   (e.g. "tumor_cell_type"). If \code{NULL}, no faceting is applied.
#' @param gene_expr_col Optional gene-level expression column name (e.g. "Cdkn2c").
#'   If not \code{NULL}, only cells with non-missing values are kept.
#' @param assay Assay from which to extract isoform data. Default: "ISO".
#' @param slot Slot to use ("data", "counts", etc.). Default: "data".
#' @param filter_detected Logical; if TRUE, keep only cells where all selected
#'   isoforms have expression > 0. Default: TRUE.
#' @param add_pvalues Logical; if TRUE, add p-values (Wilcoxon test) using
#'   \pkg{ggpubr::stat_compare_means()} if \pkg{ggpubr} >= 0.6.2 is available.
#'   Default: TRUE.
#'
#' @return A \code{ggplot} object with grouped boxplots, optional facets,
#'   and optional significance annotations.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(Matrix)
#'
#' set.seed(42)
#'
#' # --- simulate Seurat object ---
#' iso_names <- c("GeneX..Tx1", "GeneX..Tx2")
#' mat <- Matrix(
#'   rpois(200, lambda = 10),
#'   nrow = 2,
#'   sparse = TRUE,
#'   dimnames = list(iso_names, paste0("Cell", 1:100))
#' )
#'
#' obj <- CreateSeuratObject(counts = mat, assay = "ISO")
#'
#' # --- add metadata ---
#' obj$Phase <- sample(c("G1", "S", "G2M"), 100, TRUE)
#' obj$tumor_cell_type <- sample(c("Schwann", "Fibroblast", "Immune"), 100, TRUE)
#'
#' # --- test function ---
#' IsoformBoxPlot(
#'   seurat_obj  = obj,
#'   isoforms    = iso_names,
#'   group_col   = "Phase",
#'   facet_col   = "tumor_cell_type",
#'   assay       = "ISO",
#'   slot        = "counts",
#'   filter_detected = TRUE,
#'   add_pvalues = TRUE
#' )
#' }
#'
#' @importFrom Seurat GetAssayData FetchData
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot facet_wrap labs theme_bw
#' @importFrom ggplot2 theme element_text element_rect element_blank
#' @export
IsoformBoxPlot <- function(
    seurat_obj,
    isoforms,
    group_col,
    facet_col = NULL,
    gene_expr_col = NULL,
    assay = "ISO",
    slot = "data",
    filter_detected = TRUE,
    add_pvalues = TRUE
) {
  # Basic checks
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }
  meta <- seurat_obj@meta.data
  if (!group_col %in% colnames(meta)) {
    stop("group_col '", group_col, "' not found in metadata.")
  }
  if (!is.null(facet_col) && !facet_col %in% colnames(meta)) {
    stop("facet_col '", facet_col, "' not found in metadata.")
  }
  if (!all(isoforms %in% rownames(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)))) {
    missing_iso <- setdiff(
      isoforms,
      rownames(Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot))
    )
    stop("Some isoforms not found in assay '", assay, "': ",
         paste(missing_iso, collapse = ", "))
  }

  # Extract isoform matrix
  mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)[isoforms, , drop = FALSE]

  # Build data.frame with isoforms + metadata
  df <- as.data.frame(t(as.matrix(mat)))
  df[[group_col]] <- meta[[group_col]]

  if (!is.null(facet_col)) {
    df[[facet_col]] <- meta[[facet_col]]
  }

  if (!is.null(gene_expr_col)) {
    gene_vals <- Seurat::FetchData(seurat_obj, vars = gene_expr_col, slot = slot)[, 1]
    df[[gene_expr_col]] <- gene_vals
  }

  # Filter: all isoforms detected > 0
  if (filter_detected) {
    keep <- apply(df[, isoforms, drop = FALSE], 1L, function(x) all(x > 0))
    df <- df[keep, , drop = FALSE]
  }

  # Long format
  df_long <- tidyr::pivot_longer(
    df,
    cols = isoforms,
    names_to = "Isoform",
    values_to = "Expression"
  )

  # Drop rows with missing grouping
  df_long <- df_long[!is.na(df_long[[group_col]]), , drop = FALSE]

  # Base plot
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = .data[[group_col]], y = Expression, fill = Isoform)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    ggplot2::labs(
      x = group_col,
      y = "Expression",
      fill = "Isoform"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Facet if requested
  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free_y")
  }

  # Add p-values per facet/group (Wilcoxon between isoforms)
  if (add_pvalues) {
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      warning("Package 'ggpubr' not installed: p-values will not be added.")
    } else if (utils::packageVersion("ggpubr") < "0.6.2") {
      warning("ggpubr < 0.6.2 detected: skipping stat_compare_means(). Please update ggpubr to >= 0.6.2.")
    } else {
      p <- p +
        ggpubr::stat_compare_means(
          mapping = ggplot2::aes(group = Isoform),
          method  = "wilcox.test",
          label   = "p.signif",
          hide.ns = TRUE
        )
    }
  }

  p
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("Expression", "Isoform"))
}
