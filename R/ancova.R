#' Variance of the ANCOVA estimator with run-in
#' mean as covariate
#'
#' Computes the per-pair variance of the treatment
#' effect estimator from a model that regresses the
#' post-randomization mean on treatment group and the
#' run-in mean:
#'   Y_post_bar = alpha + gamma * g_i + phi * Y_pre_bar + e_i
#'
#' The variance is Var(gamma) = Var(Y_post_bar) *
#' (1 - rho^2) * 2/m, where rho is the correlation
#' between Y_pre_bar and Y_post_bar under the marginal
#' model.
#'
#' When r = 0 there is no run-in covariate and the
#' ANCOVA reduces to the simple t-test on post means.
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close
#'   observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric scalar: per-pair Var(gamma_hat).
#' @export
var_gamma_ancova <- function(r, k, f, sigma2,
                             sigma_b2, t_interval) {
  p <- r + k + f
  R <- build_R_matrix(p, sigma2)
  des <- build_design(r, k, f, t_interval)
  Z <- des$Z
  G <- matrix(sigma_b2, 1, 1)
  Sigma <- R + Z %*% G %*% t(Z)

  n_post <- k + f

  if (r == 0) {
    w_post <- rep(1 / n_post, p)
    var_post <- as.numeric(t(w_post) %*% Sigma %*% w_post)
    2 * var_post
  } else {
    w_pre <- c(rep(1 / r, r), rep(0, k + f))
    w_post <- c(rep(0, r), rep(1 / n_post, k + f))

    var_pre <- as.numeric(t(w_pre) %*% Sigma %*% w_pre)
    var_post <- as.numeric(t(w_post) %*% Sigma %*% w_post)
    cov_pre_post <- as.numeric(
      t(w_pre) %*% Sigma %*% w_post)

    rho_sq <- cov_pre_post^2 / (var_pre * var_post)
    2 * var_post * (1 - rho_sq)
  }
}


#' Relative efficiency of ANCOVA vs GLS estimator
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Post-randomization observations.
#' @param f Integer. Common close observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric: RE = Var(GLS) / Var(ANCOVA).
#' @export
re_ancova_vs_gls <- function(r, k, f, sigma2,
                             sigma_b2, t_interval) {
  v_gls <- var_gamma_matrix(r, k, f, sigma2, sigma_b2,
                            t_interval)
  v_anc <- var_gamma_ancova(r, k, f, sigma2, sigma_b2,
                            t_interval)
  v_gls / v_anc
}
