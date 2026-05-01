#' Marginal variance reduction from the r-th run-in
#' observation
#'
#' Computes V(r-1, k, f) - V(r, k, f), the reduction
#' in per-pair variance from adding one additional
#' run-in observation.
#'
#' @param r Integer. Run-in observation being added
#'   (must be >= 1).
#' @param k Integer. Post-randomization observations.
#' @param f Integer. Common close observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Numeric scalar: variance reduction.
#' @export
marginal_reduction <- function(r, k, f, sigma2,
                               sigma_b2, t_interval) {
  stopifnot(r >= 1)
  var_gamma_matrix(r - 1, k, f, sigma2, sigma_b2,
                   t_interval) -
    var_gamma_matrix(r, k, f, sigma2, sigma_b2,
                     t_interval)
}


#' Optimal allocation of observations across phases
#'
#' Given a fixed total number of change scores p,
#' finds the allocation (r, k, f) that minimizes
#' Var(gamma_hat) subject to r + k + f = p, k >= 1.
#'
#' @param p Integer. Total number of change scores.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return A list with components r, k, f, variance.
#' @export
optimal_allocation <- function(p, sigma2, sigma_b2,
                               t_interval) {
  stopifnot(p >= 2)
  configs <- expand.grid(
    r = 0:(p - 1),
    k = 1:p,
    f = 0:(p - 1)
  )
  configs <- configs[configs$r + configs$k +
                       configs$f == p, ]

  configs$variance <- sapply(seq_len(nrow(configs)),
    function(i) {
      var_gamma_matrix(
        configs$r[i], configs$k[i], configs$f[i],
        sigma2, sigma_b2, t_interval
      )
    })

  best <- configs[which.min(configs$variance), ]
  list(
    r = best$r,
    k = best$k,
    f = best$f,
    variance = best$variance
  )
}


#' Common close benefit threshold
#'
#' For a given (r, k), determines whether adding one
#' common close observation reduces the variance.
#'
#' @param r Integer. Run-in observations.
#' @param k Integer. Post-randomization observations.
#' @param sigma2 Numeric. Residual variance.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @return Logical: TRUE if f=1 improves over f=0.
#' @export
common_close_beneficial <- function(r, k, sigma2,
                                    sigma_b2,
                                    t_interval) {
  v0 <- var_gamma_matrix(r, k, 0, sigma2, sigma_b2,
                         t_interval)
  v1 <- var_gamma_matrix(r, k, 1, sigma2, sigma_b2,
                         t_interval)
  v1 < v0
}
