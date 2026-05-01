#' Build design matrices for a run-in clinical trial
#'
#' Constructs the fixed effects design matrices (treatment
#' and placebo) and the random effects design vector for
#' a trial with r run-in, k treatment, and f common close
#' observations.
#'
#' When r = 0, the model has two fixed effects (beta,
#' gamma). When r > 0, a third parameter delta (run-in
#' slope) is added. The gamma column encodes the
#' accumulated treatment effect: linear in the treatment
#' period, constant (at k) during the common close.
#'
#' @param r Integer. Number of run-in observations.
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close observations.
#' @param t_interval Numeric. Observation spacing.
#' @return A list with components:
#'   \describe{
#'     \item{X_trt}{Fixed effects design matrix for
#'       treatment.}
#'     \item{X_plc}{Fixed effects design matrix for
#'       placebo.}
#'     \item{Z}{Random effects design vector (p x 1).}
#'     \item{p}{Total number of change scores.}
#'     \item{ncol_X}{Number of fixed effects.}
#'   }
#' @export
build_design <- function(r, k, f, t_interval) {
  p <- r + k + f
  stopifnot(p >= 2)

  d_runin <- if (r > 0) seq(-r, -1) else numeric(0)
  d_trt <- if (k > 0) seq(1, k) else numeric(0)
  d_close <- if (f > 0) seq(k + 1, k + f) else numeric(0)
  d_all <- c(d_runin, d_trt, d_close)

  Z <- t_interval * d_all

  gamma_trt <- c(rep(0, r), seq(1, k), rep(k, f))
  gamma_plc <- rep(0, p)

  if (r == 0) {
    X_trt <- cbind(
      beta  = t_interval * c(d_trt, d_close),
      gamma = t_interval * gamma_trt
    )
    X_plc <- cbind(
      beta  = t_interval * c(d_trt, d_close),
      gamma = t_interval * gamma_plc
    )
  } else {
    X_trt <- cbind(
      delta = t_interval * c(d_runin, rep(0, k + f)),
      beta  = t_interval * c(rep(0, r), d_trt, d_close),
      gamma = t_interval * gamma_trt
    )
    X_plc <- cbind(
      delta = t_interval * c(d_runin, rep(0, k + f)),
      beta  = t_interval * c(rep(0, r), d_trt, d_close),
      gamma = t_interval * gamma_plc
    )
  }

  list(
    X_trt = X_trt,
    X_plc = X_plc,
    Z = matrix(Z, ncol = 1),
    p = p,
    ncol_X = ncol(X_trt)
  )
}
