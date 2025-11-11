#' Plot isoform PI or expression along pseudotime with annotation bands
#'
#' Visualizes isoform-level usage or expression for a given gene along pseudotime,
#' together with one or multiple categorical annotations displayed as horizontal
#' bands (e.g. cell type, cell cycle phase). The two most abundant isoforms of
#' the gene are shown.
#'
#' @param seurat_obj A Seurat object.
#' @param gene Character scalar. Gene name used as prefix for isoform IDs
#'   (e.g. "GeneA" for rows like "GeneA..Tx1").
#' @param pseudotime_col Name of the metadata column containing pseudotime values.
#' @param band_cols Character vector of metadata column names to display as
#'   annotation bands (e.g. \code{c("cell_type", "Phase")}). Each should be
#'   a categorical variable.
#' @param assay Assay name containing isoform-level counts. Default: \code{"ISO"}.
#' @param slot Slot of the assay to use (e.g. \code{"counts"}, \code{"data"}).
#'   Default: \code{"counts"}.
#' @param mode Either \code{"PI"} (default) to plot isoform Percent Inclusion
#'   (per-cell proportions across isoforms of the gene) or \code{"exp"} to plot
#'   raw isoform expression.
#' @param delimiter Regex delimiter between gene and isoform ID.
#'   Default: \code{"\\\\.\\\\."}.
#' @param nbins Integer. Number of bins along pseudotime for aggregating
#'   annotation bands. Default: 100.
#'
#' @return A \code{patchwork} object combining the annotation band plot(s)
#'   and the isoform trend plot.
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(Matrix)
#'
#' set.seed(123)
#'
#' iso_names <- c("GeneA..Tx1", "GeneA..Tx2", "GeneA..Tx3")
#' mat <- Matrix(
#'   rpois(300, lambda = 10),
#'   nrow = 3,
#'   sparse = TRUE,
#'   dimnames = list(iso_names, paste0("Cell", 1:100))
#' )
#'
#' obj <- CreateSeuratObject(counts = mat, assay = "ISO")
#' obj$pseudotime <- sort(runif(100, 0, 1))
#' obj$cell_type  <- sample(c("Tumor", "Stromal", "Immune"), 100, TRUE)
#' obj$Phase      <- sample(c("G1", "S", "G2M"), 100, TRUE)
#'
#' IsoformPseudotimePlot(
#'   seurat_obj     = obj,
#'   gene           = "GeneA",
#'   pseudotime_col = "pseudotime",
#'   band_cols      = c("cell_type", "Phase"),
#'   assay          = "ISO",
#'   mode           = "PI"
#' )
#' }
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr %>% group_by summarise arrange mutate filter all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_tile geom_smooth
#' @importFrom ggplot2 scale_color_manual scale_fill_manual
#' @importFrom ggplot2 labs theme_minimal theme theme_void guides guide_legend
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stringr str_replace
#' @importFrom scales hue_pal
#' @importFrom stats setNames
#' @export
IsoformPseudotimePlot <- function(
    seurat_obj,
    gene,
    pseudotime_col,
    band_cols = "cell_type",
    assay = "ISO",
    slot = "counts",
    mode = c("PI", "exp"),
    delimiter = "\\.\\.",
    nbins = 100
) {
  mode <- match.arg(mode)

  # ---- Checks ----
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }

  meta <- seurat_obj@meta.data

  if (!pseudotime_col %in% colnames(meta)) {
    stop("Column '", pseudotime_col, "' not found in seurat_obj@meta.data.")
  }
  if (!all(band_cols %in% colnames(meta))) {
    missing_cols <- setdiff(band_cols, colnames(meta))
    stop("band_cols not found in metadata: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(nbins) || length(nbins) != 1L || nbins <= 0) {
    stop("nbins must be a single positive numeric value.")
  }

  # ---- 1. Isoforms of the gene ----
  assay_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  all_features <- rownames(assay_mat)
  pattern <- paste0("^", gene, delimiter)
  genes <- all_features[grepl(pattern, all_features)]

  if (length(genes) < 2L) {
    stop("Less than 2 isoforms found for gene: ", gene)
  }

  expr_data <- assay_mat[genes, , drop = FALSE]

  pseudotime <- meta[[pseudotime_col]]
  if (all(is.na(pseudotime))) {
    stop("All pseudotime values are NA in column '", pseudotime_col, "'.")
  }

  # ---- 2. Compute PI or expression ----
  mat <- as.matrix(expr_data)  # isoforms x cells

  if (mode == "exp") {
    value_mat <- mat
    y_lab <- "Expression"
  } else {
    gene_totals <- base::colSums(mat)
    gene_totals[gene_totals == 0] <- NA_real_
    value_mat <- sweep(mat, 2, gene_totals, "/")
    value_mat[!is.finite(value_mat)] <- NA_real_
    y_lab <- "PI"
  }

  # ---- 3. Build long format & select top2 isoforms ----
  df_val <- data.frame(
    pseudotime = pseudotime,
    meta[, band_cols, drop = FALSE],
    t(value_mat),
    check.names = FALSE
  )
  df_val <- df_val[is.finite(df_val$pseudotime), , drop = FALSE]

  df_long <- tidyr::pivot_longer(
    df_val,
    cols = dplyr::all_of(genes),
    names_to = "isoform",
    values_to = "value"
  )

  df_long$isoform_simple <- stringr::str_replace(
    df_long$isoform,
    paste0("^", gene, delimiter),
    ""
  )

  iso_means <- df_long %>%
    dplyr::group_by(isoform_simple) %>%
    dplyr::summarise(
      mean_val = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_val))

  top2 <- iso_means$isoform_simple[seq_len(min(2L, nrow(iso_means)))]
  if (length(top2) == 0L) {
    stop("No isoforms with valid values to plot.")
  }

  pal_iso <- stats::setNames(
    c("#E41A1C", "#377EB8")[seq_along(top2)],
    top2
  )

  df_long <- df_long %>%
    dplyr::filter(isoform_simple %in% top2)

  # ---- 4. Helper: one band per annotation column ----
  # Palette auto différente par bande : chaque bande utilise une gamme de teintes décalée.
  make_band_plot <- function(pseudotime_vec, annotation_vec, band_name, nbins, band_index) {
    df_band <- data.frame(
      pseudotime = pseudotime_vec,
      value = annotation_vec
    )
    df_band <- df_band[
      is.finite(df_band$pseudotime) & !is.na(df_band$value),
      ,
      drop = FALSE
    ]
    if (nrow(df_band) == 0L) {
      return(NULL)
    }

    df_band$pt_bin <- cut(
      df_band$pseudotime,
      breaks = nbins,
      include.lowest = TRUE
    )

    df_band_binned <- df_band %>%
      dplyr::group_by(pt_bin) %>%
      dplyr::summarise(
        pseudotime = mean(pseudotime),
        value = {
          tab <- table(value)
          names(tab)[which.max(tab)]
        },
        .groups = "drop"
      )

    vals <- unique(df_band_binned$value)
    ntypes <- length(vals)
    if (ntypes == 0L) return(NULL)

    # Palette spécifique à cette bande (shift de teinte selon band_index)
    # On décale la gamme de teintes pour éviter de réutiliser exactement les mêmes couleurs
    hue_start <- 15 + (band_index - 1L) * 90
    hue_end   <- hue_start + 270
    pal_vals <- scales::hue_pal(h = c(hue_start, hue_end))(ntypes)
    names(pal_vals) <- vals

    ggplot2::ggplot(
      df_band_binned,
      ggplot2::aes(x = pseudotime, y = 1, fill = value)
    ) +
      ggplot2::geom_tile(
        width = diff(range(df_band$pseudotime, na.rm = TRUE)) / nbins,
        height = 1
      ) +
      ggplot2::scale_fill_manual(values = pal_vals) +
      ggplot2::theme_void() +
      ggplot2::guides(
        fill = ggplot2::guide_legend(title = band_name)
      )
  }

  meta_band <- data.frame(
    pseudotime = pseudotime,
    meta[, band_cols, drop = FALSE],
    check.names = FALSE
  )

  band_plots <- lapply(seq_along(band_cols), function(i) {
    bc <- band_cols[i]
    make_band_plot(
      pseudotime_vec = meta_band$pseudotime,
      annotation_vec = meta_band[[bc]],
      band_name = bc,
      nbins = nbins,
      band_index = i
    )
  })
  band_plots <- band_plots[!vapply(band_plots, is.null, logical(1L))]

  # ---- 5. Isoform trend plot ----
  p_iso <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = pseudotime,
      y = value,
      color = isoform_simple,
      group = isoform_simple
    )
  ) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, size = 1.1) +
    ggplot2::scale_color_manual(values = pal_iso) +
    ggplot2::labs(
      x = "Pseudotime",
      y = y_lab,
      color = "Isoform",
      title = paste0(gene, " - top isoforms: ", paste(top2, collapse = ", "))
    ) +
    ggplot2::theme_minimal()

  # ---- 6. Combine bands + main plot ----
  plots <- c(band_plots, list(p_iso))
  heights <- c(rep(0.06, length(band_plots)), 1)

  patchwork::wrap_plots(
    plots,
    ncol = 1,
    heights = heights
  ) +
    patchwork::plot_layout(guides = "collect")
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "isoform_simple",
    "mean_val",
    "pt_bin",
    "value"
  ))
}
