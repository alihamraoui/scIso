#' Plot mean PI vs Phi for isoforms
#'
#' This function summarizes isoform-level PI values and visualizes the relationship
#' between mean PI and Phi across isoforms. Isoforms can be filtered based on
#' their mean PI and total gene-level support before plotting.
#'
#' @param pi_mat A numeric or sparse matrix of PI values (isoforms × cells).
#' @param iso_mat A numeric or sparse matrix of isoform-level counts
#'   (same dimensions and rownames as \code{pi_mat}), used to derive gene-level
#'   coverage.
#' @param pi_min Minimum mean PI threshold to retain isoforms (default: 0.3).
#' @param pi_max Maximum mean PI threshold to retain isoforms (default: 0.7).
#' @param min_gene_count Minimum total gene count required to include an isoform
#'   in the plot (default: 3500).
#'
#' @details
#' For each isoform, the function computes:
#' \itemize{
#'   \item mean PI across cells;
#'   \item Phi using \code{\link{computePHI}} (inter-cell variability);
#'   \item number of cells with non-zero counts;
#'   \item total gene-level count (sum across all isoforms of a gene).
#' }
#'
#' Isoforms are filtered according to \code{pi_min}, \code{pi_max}, and
#' \code{min_gene_count}, and only those passing all filters are plotted.
#'
#' @return A \code{ggplot} object visualizing mean PI vs Phi, colored by
#' gene-level count.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' pi_mat <- matrix(runif(200, 0.05, 0.95), nrow = 50,
#'                  dimnames = list(paste0("Gene", rep(1:10, each = 5), "..Tx", 1:5),
#'                                  paste0("Cell", 1:4)))
#' iso_mat <- matrix(sample(0:20, 200, replace = TRUE), nrow = 50,
#'                   dimnames = dimnames(pi_mat))
#'
#' # Default thresholds
#' plotPIPHI(pi_mat, iso_mat)
#'
#' # Custom thresholds
#' plotPIPHI(pi_mat, iso_mat, pi_min = 0.2, pi_max = 0.8, min_gene_count = 10)
#' }
#'
#' @export
plotPIPHI <- function(pi_mat, iso_mat,
                      pi_min = 0.3,
                      pi_max = 0.7,
                      min_gene_count = 1000) {
  if (is.null(dim(pi_mat))) {
    stop("pi_mat must be a matrix.")
  }
  if (is.null(dim(iso_mat))) {
    stop("iso_mat must be a matrix.")
  }
  if (!identical(base::rownames(pi_mat), base::rownames(iso_mat))) {
    stop("pi_mat and iso_mat must have identical rownames (isoforms).")
  }
  if (!identical(base::colnames(pi_mat), base::colnames(iso_mat))) {
    stop("pi_mat and iso_mat must have identical colnames (cells).")
  }

  phi_vec <- computePHI(pi_mat)

  mean_PI <- base::data.frame(
    Isoform   = base::rownames(pi_mat),
    mean_pi   = base::rowMeans(pi_mat, na.rm = TRUE),
    phi       = phi_vec,
    count     = Matrix::rowSums(iso_mat != 0),
    stringsAsFactors = FALSE
  )

  mean_PI$gene <- stringr::str_extract(mean_PI$Isoform, "^[^\\.]+")

  gene_counts <- mean_PI |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      gene_count = base::sum(count, na.rm = TRUE),
      .groups = "drop"
    )

  mean_PI <- dplyr::left_join(mean_PI, gene_counts, by = "gene")

  data_mid <- mean_PI[
    mean_PI$mean_pi > pi_min &
      mean_PI$mean_pi < pi_max &
      mean_PI$gene_count > min_gene_count &
      !base::is.na(mean_PI$phi),
    ,
    drop = FALSE
  ]

  if (nrow(data_mid) == 0L) {
    base::warning(data_mid)#"No isoforms passed the PI/Phi filtering criteria.")
    return(ggplot2::ggplot())
  }

  data_mid$Gene <- base::sub("\\.\\..*$", "", data_mid$Isoform)

  ggplot2::ggplot(
    data_mid,
    ggplot2::aes(x = mean_pi, y = phi, color = gene_count)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_gradientn(
      colors = base::rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$Gene),
      size = 3
    ) +
    ggplot2::xlab(latex2exp::TeX("$\\bar{\\pi}$")) +
    ggplot2::ylab(latex2exp::TeX("$\\phi$")) +
    ggplot2::labs(color = "Gene count") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 15)
    )
}
