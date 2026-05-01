#' Variance of the treatment effect estimator via
#' Woodbury identity
#'
#' Computes Var(gamma_hat) for a single participant
#' pair (one treatment, one placebo) using the matrix
#' form of the Woodbury identity:
#' X' Sigma^{-1} X = X' R^{-1} X - W, where W is the
#' sum of per-group rank-one Woodbury corrections.
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close observations.
#' @param sigma2 Numeric. Residual variance (sigma^2).
#' @param sigma_b2 Numeric. Random slope variance
#'   (sigma_b^2).
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric scalar: per-pair Var(gamma_hat).
#' @export
var_gamma_matrix <- function(r, k, f, sigma2, sigma_b2,
                             t_interval) {
  p <- r + k + f
  R <- build_R_matrix(p, sigma2)
  Ri <- solve(R)
  des <- build_design(r, k, f, t_interval)
  Z <- des$Z
  G <- matrix(sigma_b2, 1, 1)

  H <- solve(solve(G) + t(Z) %*% Ri %*% Z)

  Xt <- des$X_trt
  Xp <- des$X_plc
  XRX <- t(Xt) %*% Ri %*% Xt + t(Xp) %*% Ri %*% Xp

  XtRZ <- t(Xt) %*% Ri %*% Z
  XpRZ <- t(Xp) %*% Ri %*% Z
  W <- XtRZ %*% H %*% t(XtRZ) + XpRZ %*% H %*% t(XpRZ)

  XSiX <- XRX - W
  V <- solve(XSiX)
  nc <- des$ncol_X
  V[nc, nc]
}


#' Closed-form variance for the standard design
#' (Frost et al. 2008, Eq. 13)
#'
#' Per-pair variance with r = 0, k = 2, f = 0.
#'
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric scalar.
#' @export
var_gamma_frost <- function(sigma2, sigma_b2, t_interval) {
  (sigma2 + 2 * sigma_b2 * t_interval^2) / t_interval^2
}


#' Closed-form variance for the two-run-in design
#' (r = 2, k = 1, f = 0)
#'
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric scalar.
#' @export
var_gamma_r2 <- function(sigma2, sigma_b2, t_interval) {
  D <- sigma2 + 5 * sigma_b2 * t_interval^2
  8 * sigma2 * D /
    (3 * t_interval^2 *
       (sigma2 + 2 * sigma_b2 * t_interval^2))
}
