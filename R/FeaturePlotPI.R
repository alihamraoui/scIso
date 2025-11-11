#' Visualize isoform percent inclusion (PI), ratios and distributions on embeddings
#'
#' This function visualizes the relative usage (PI, Percent Inclusion) of two
#' isoforms across cells on a low-dimensional embedding stored in a Seurat object.
#' It computes per-cell PI for each isoform, their ratio and log2-ratio, adds
#' these values to the object metadata, and displays:
#'
#' \itemize{
#'   \item PI of isoform 1 on the embedding,
#'   \item PI of isoform 2 on the embedding,
#'   \item log2(ratio iso1 / iso2) on the embedding,
#'   \item histogram of the PI distribution for isoform 1.
#' }
#'
#' @import patchwork
#'
#' @param obj A Seurat object.
#' @param iso1 Isoform ID for the first isoform (must be a row name in the assay).
#' @param iso2 Isoform ID for the second isoform (must be a row name in the assay).
#' @param assay Assay name containing isoform counts. Default: \code{"ISO"}.
#' @param slot Slot to use in the assay (e.g. \code{"counts"}, \code{"data"}).
#'   Default: \code{"counts"}.
#' @param reduction Name of the dimensional reduction to use (e.g. \code{"umap"}).
#' @param name Optional base name used to construct metadata column names.
#' @param add_ratio Logical; if \code{TRUE}, add ratio and log2-ratio to metadata.
#'   Default: \code{TRUE}.
#' @param min_total Numeric threshold: cells with \code{iso1 + iso2 < min_total}
#'   are set to \code{NA} for PI. Default: \code{0}.
#' @param col Suffix used for PI-based metadata/label naming. Default: \code{"PI"}.
#' @param bins Number of bins for the PI histogram. Default: \code{30}.
#'
#' @return A \code{patchwork} object combining three feature plots and one histogram.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(Matrix)
#'
#' set.seed(42)
#'
#' mat <- Matrix(
#'   rpois(200, lambda = 10),
#'   nrow = 2,
#'   sparse = TRUE,
#'   dimnames = list(
#'     c("GeneA..Tx1", "GeneA..Tx2"),
#'     paste0("Cell", 1:100)
#'   )
#' )
#'
#' obj <- CreateSeuratObject(counts = mat, assay = "ISO")
#'
#' emb <- matrix(
#'   rnorm(200),
#'   ncol = 2,
#'   dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2"))
#' )
#'
#' obj[["umap"]] <- CreateDimReducObject(
#'   embeddings = emb,
#'   key = "UMAP_",
#'   assay = "ISO"
#' )
#'
#' FeaturePlotPI(obj, "GeneA..Tx1", "GeneA..Tx2", reduction = "umap")
#' }
#'
#' @export
FeaturePlotPI <- function(
    obj,
    iso1,
    iso2,
    reduction,
    assay = "ISO",
    slot  = "counts",
    name  = NULL,
    add_ratio = TRUE,
    min_total = 0,
    col = "PI",
    bins = 30
){
  stopifnot(inherits(obj, "Seurat"))
  if (!assay %in% names(obj@assays)) stop("Assay introuvable: ", assay)

  mat <- Seurat::GetAssayData(obj, assay = assay, slot = slot)
  miss <- base::setdiff(c(iso1, iso2), base::rownames(mat))
  if (length(miss) > 0) {
    stop(
      "Isoforme(s) absente(s) de l'assay '", assay, "': ",
      paste(miss, collapse = ", ")
    )
  }

  v1 <- base::as.numeric(mat[iso1, ])
  v2 <- base::as.numeric(mat[iso2, ])
  tot <- v1 + v2

  pi1 <- ifelse(tot > 0, v1 / tot, NA_real_)
  pi1[tot < min_total] <- NA_real_

  pi2 <- ifelse(tot > 0, v2 / tot, NA_real_)
  pi2[tot < min_total] <- NA_real_

  ratio     <- ifelse(v2 > 0, v1 / v2, NA_real_)
  log2ratio <- base::log2((v1 + 1) / (v2 + 1))
  log2ratio[log2ratio == 0] <- NA_real_

  base_name <- if (!is.null(name)) {
    name
  } else {
    base::gsub("[^A-Za-z0-9]+", "_", paste0("ISO_", iso1, "_vs_", iso2))
  }

  trans1 <- stringr::str_split(iso1, "\\.\\.")[[1]][2]
  trans2 <- stringr::str_split(iso2, "\\.\\.")[[1]][2]

  obj[[paste0(trans1, "_PI")]] <- pi1
  obj[[paste0(trans2, "_PI")]] <- pi2

  if (add_ratio) {
    obj[[paste0(base_name, "_ratio")]]     <- ratio
    obj[[paste0(base_name, "_log2ratio")]] <- log2ratio
  }

  # Palettes (utilisées via `cols =` dans FeaturePlot)
  pi_pal   <- paletteer::paletteer_c("ggthemes::Orange-Blue-White Diverging", 30)
  ratio_pal <- RColorBrewer::brewer.pal(11, "Spectral")

  # 1) PI isoform 1
  p_iso1 <- Seurat::FeaturePlot(
    obj,
    reduction = reduction,
    features = paste0(trans1, "_PI"),
    cols = base::rev(pi_pal)
  ) +
    ggplot2::labs(color = "PI") +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    Seurat::NoAxes()

  # 2) PI isoform 2
  p_iso2 <- Seurat::FeaturePlot(
    obj,
    reduction = reduction,
    features = paste0(trans2, "_PI"),
    cols = base::rev(pi_pal)
  ) +
    ggplot2::labs(color = "PI") +
    ggplot2::theme(
      aspect.ratio = 1,
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    Seurat::NoAxes()

  # 3) log2 ratio
  p_log2ratio <- Seurat::FeaturePlot(
    obj,
    reduction = reduction,
    features = paste0(base_name, "_log2ratio"),
    cols = ratio_pal
  ) +
    ggplot2::labs(title = "log2ratio", color = "log2ratio")+
    ggplot2::theme(
      aspect.ratio = 1,
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    Seurat::NoAxes()

  # 4) Histogram of PI for iso1
  df <- obj@meta.data

  p_hist <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[paste0(trans1, "_PI")]])
  ) +
    ggplot2::geom_histogram(bins = bins, na.rm = TRUE) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = col,
      y = "Number of cells",
      title = "PI across cells"
    )

  p_iso1 | p_iso2 | p_log2ratio | p_hist
}
