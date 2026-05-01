#' Required sample size per group
#'
#' Computes the number of participants per group for a
#' two-sided test at level alpha with the specified power,
#' based on the per-pair variance from
#' [var_gamma_matrix()].
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @param delta Numeric. Minimum clinically important
#'   treatment effect on the rate of change.
#' @param alpha Numeric. Two-sided significance level
#'   (default 0.05).
#' @param power Numeric. Desired power (default 0.90).
#' @return Integer: required sample size per group.
#' @export
sample_size <- function(r, k, f, sigma2, sigma_b2,
                        t_interval, delta, alpha = 0.05,
                        power = 0.90) {
  v1 <- var_gamma_matrix(r, k, f, sigma2, sigma_b2,
                         t_interval)
  za <- qnorm(1 - alpha / 2)
  zb <- qnorm(power)
  ceiling((za + zb)^2 * v1 / delta^2)
}


#' Statistical power for a given sample size
#'
#' Computes the power of a two-sided test at level alpha
#' for a given number of participants per group.
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @param delta Numeric. Treatment effect.
#' @param m Integer. Number of participants per group.
#' @param alpha Numeric. Two-sided significance level
#'   (default 0.05).
#' @return Numeric: statistical power.
#' @export
power_calc <- function(r, k, f, sigma2, sigma_b2,
                       t_interval, delta, m,
                       alpha = 0.05) {
  v1 <- var_gamma_matrix(r, k, f, sigma2, sigma_b2,
                         t_interval)
  za <- qnorm(1 - alpha / 2)
  pnorm(delta * sqrt(m) / sqrt(v1) - za)
}


#' Relative efficiency of a run-in design
#'
#' Ratio of the standard design variance to the run-in
#' design variance. Values greater than 1 indicate the
#' run-in design requires fewer participants.
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric: relative efficiency.
#' @export
relative_efficiency <- function(r, k, f, sigma2,
                                sigma_b2, t_interval) {
  v_base <- var_gamma_matrix(0, k, 0, sigma2, sigma_b2,
                             t_interval)
  v_design <- var_gamma_matrix(r, k, f, sigma2, sigma_b2,
                               t_interval)
  v_base / v_design
}
