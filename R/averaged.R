#' Variance of the averaged change-score estimator
#'
#' Computes the per-pair variance of the simple
#' difference-of-means estimator that averages all
#' post-randomization change scores and subtracts the
#' average of all run-in change scores. This estimator
#' does not use the covariance structure and is
#' analyzed via a t-test.
#'
#' Under the model, the averaged change score for
#' participant i is:
#'   d_i = (1/c_i) sum(post) - (1/r_i) sum(pre)
#' where c_i = k + f and r_i = r.
#'
#' The variance is computed from the marginal
#' covariance Sigma = R + ZGZ' by applying the
#' contrast vector that computes d_i.
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close
#'   observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric scalar: per-pair Var(d_trt - d_plc).
#' @export
var_gamma_avg <- function(r, k, f, sigma2, sigma_b2,
                          t_interval) {
  p <- r + k + f
  R <- build_R_matrix(p, sigma2)
  des <- build_design(r, k, f, t_interval)
  Z <- des$Z
  G <- matrix(sigma_b2, 1, 1)
  Sigma <- R + Z %*% G %*% t(Z)

  n_post <- k + f
  if (r == 0) {
    w_trt <- rep(1 / n_post, p)
    w_plc <- rep(1 / n_post, p)
  } else {
    w_pre <- c(rep(-1 / r, r), rep(0, k + f))
    w_post <- c(rep(0, r), rep(1 / n_post, k + f))
    w_trt <- w_pre + w_post
    w_plc <- w_trt
  }

  var_d <- as.numeric(t(w_trt) %*% Sigma %*% w_trt)
  2 * var_d
}


#' Relative efficiency of averaged vs GLS estimator
#'
#' Ratio of GLS variance to averaged variance.
#' Values less than 1 indicate GLS is more efficient
#' (which is always the case under correct
#' specification).
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Post-randomization observations.
#' @param f Integer. Common close observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric: RE = Var(GLS) / Var(avg).
#' @export
re_avg_vs_gls <- function(r, k, f, sigma2, sigma_b2,
                          t_interval) {
  v_gls <- var_gamma_matrix(r, k, f, sigma2, sigma_b2,
                            t_interval)
  v_avg <- var_gamma_avg(r, k, f, sigma2, sigma_b2,
                         t_interval)
  v_gls / v_avg
}
