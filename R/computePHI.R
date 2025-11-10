#' Compute Phi (ϕ) from a PI matrix
#'
#' This function estimates the dispersion parameter \eqn{\phi} of a Beta
#' distribution for each isoform (row) in a matrix of PI (Percent Inclusion)
#' values. The parameter \eqn{\phi} reflects the intercellular variability of
#' isoform usage across cells.
#'
#' @param pi_mat A numeric or sparse matrix of PI values (isoforms × cells).
#'
#' @details
#' For each isoform, \eqn{\phi} is estimated from the sample mean (\eqn{m}) and
#' variance (\eqn{v}) of PI values across cells, using the method of moments:
#' \deqn{\phi = 1 / (\alpha + \beta + 1)}
#' where \eqn{\alpha = m((1 - m)/v - 1)} and \eqn{\beta = (1 - m)((1 - m)/v - 1)}.
#'
#' Isoforms with fewer than three valid observations or with degenerate mean or
#' variance values return \code{NA}.
#'
#' @return A named numeric vector of \eqn{\phi} values, one per isoform (row of
#' \code{pi_mat}). Rows with insufficient data or invalid parameters return
#' \code{NA}.
#'
#' @examples
#' set.seed(123)
#' pi_mat <- matrix(runif(20, 0.1, 0.9), nrow = 5,
#'                  dimnames = list(paste0("Iso", 1:5), paste0("Cell", 1:4)))
#' phi_vec <- computePHI(pi_mat)
#' phi_vec
#'
#' @export
computePHI <- function(pi_mat) {

  # Internal function: compute phi for a single isoform vector
  compute_phi_single <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    m <- mean(x)
    v <- stats::var(x)

    # Degenerate or incompatible cases
    if (n < 3 || v == 0 || m == 0 || m == 1 || is.na(m) || is.na(v)) {
      return(NA_real_)
    }

    common <- m * (1 - m) / v - 1
    if (common <= 0) return(NA_real_)

    alpha <- m * common
    beta  <- (1 - m) * common
    if (alpha <= 0 || beta <= 0) return(NA_real_)

    1 / (alpha + beta + 1)
  }

  # Apply per isoform
  phi_values <- base::apply(pi_mat, 1, compute_phi_single)
  names(phi_values) <- base::rownames(pi_mat)

  return(phi_values)
}

