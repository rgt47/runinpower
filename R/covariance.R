#' Build residual covariance matrix for change scores
#'
#' Constructs the p x p residual covariance matrix
#' R = sigma2 * (I_p + 1_p 1_p') for change scores
#' relative to baseline. The off-diagonal structure
#' arises because all change scores share the baseline
#' error e_{i0}.
#'
#' @param p Integer. Number of change scores.
#' @param sigma2 Numeric. Residual variance.
#' @return A p x p symmetric positive definite matrix.
#' @export
build_R_matrix <- function(p, sigma2) {
  sigma2 * (diag(p) + matrix(1, p, p))
}
