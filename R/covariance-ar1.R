#' Build residual covariance matrix under AR(1)
#'
#' Constructs the p x p residual covariance matrix for
#' change scores relative to baseline when the raw
#' errors follow an AR(1) process with parameter rho.
#' Element (j, l) is Cov(e_j - e_0, e_l - e_0) =
#' sigma2 * (rho^|j-l| - rho^|j| - rho^|l| + 1).
#'
#' @param p Integer. Number of change scores.
#' @param sigma2 Numeric. Residual variance.
#' @param rho Numeric. AR(1) autocorrelation (0 < rho < 1).
#' @param d_all Integer vector of length p. Time indices
#'   relative to baseline (e.g., c(-2, -1, 1) for r=2, k=1).
#'   If NULL, defaults to change-from-baseline with equally
#'   spaced times.
#' @return A p x p symmetric positive definite matrix.
#' @export
build_R_matrix_ar1 <- function(p, sigma2, rho,
                                d_all = NULL) {
  if (is.null(d_all)) {
    stop("d_all must be provided for AR(1) covariance")
  }
  R <- matrix(0, p, p)
  for (j in seq_len(p)) {
    for (l in seq_len(p)) {
      R[j, l] <- sigma2 * (
        rho^abs(d_all[j] - d_all[l]) -
        rho^abs(d_all[j]) -
        rho^abs(d_all[l]) + 1
      )
    }
  }
  R
}


#' Variance of treatment effect under AR(1) correlation
#'
#' Computes Var(gamma_hat) per participant pair using
#' the Woodbury identity with an AR(1) residual
#' covariance structure.
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @param rho Numeric. AR(1) autocorrelation parameter.
#' @return Numeric scalar: per-pair Var(gamma_hat).
#' @export
var_gamma_ar1 <- function(r, k, f, sigma2, sigma_b2,
                          t_interval, rho) {
  p <- r + k + f
  des <- build_design(r, k, f, t_interval)

  d_runin <- if (r > 0) seq(-r, -1) else numeric(0)
  d_trt <- if (k > 0) seq(1, k) else numeric(0)
  d_close <- if (f > 0) seq(k + 1, k + f) else numeric(0)
  d_all <- c(d_runin, d_trt, d_close)

  R <- build_R_matrix_ar1(p, sigma2, rho, d_all)
  Ri <- solve(R)
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
