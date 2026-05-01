#' Variance with replicated endpoint observations
#'
#' Computes Var(gamma_hat) when the final treatment
#' visit is replaced by n_rep closely-spaced
#' replicate measurements. The replicates are assumed
#' to be at essentially the same time point (separated
#' by days, not months), so each replicate has the same
#' random slope contribution but independent
#' measurement error. The effective residual variance
#' for the averaged endpoint is sigma2 / n_rep.
#'
#' Similarly, n_rep_pre replicates can be taken at
#' the baseline visit, reducing the baseline
#' measurement error contribution.
#'
#' @param r Integer. Number of run-in observations
#'   (each at distinct time points, not replicates).
#' @param k Integer. Number of post-randomization
#'   observations.
#' @param f Integer. Number of common close
#'   observations.
#' @param sigma2 Numeric. Residual variance per
#'   single observation.
#' @param sigma_b2 Numeric. Random slope variance.
#' @param t_interval Numeric. Observation spacing.
#' @param n_rep_post Integer. Number of closely-spaced
#'   replicates at the final treatment visit
#'   (default 1, no replication).
#' @param n_rep_pre Integer. Number of closely-spaced
#'   replicates at the baseline visit
#'   (default 1, no replication).
#' @return Numeric scalar: per-pair Var(gamma_hat).
#' @export
var_gamma_replicated <- function(r, k, f, sigma2,
                                  sigma_b2, t_interval,
                                  n_rep_post = 1,
                                  n_rep_pre = 1) {
  p <- r + k + f
  des <- build_design(r, k, f, t_interval)
  Z <- des$Z
  G <- matrix(sigma_b2, 1, 1)

  sigma2_vec <- rep(sigma2, p)

  if (n_rep_post > 1 && p >= 1) {
    if (f > 0) {
      sigma2_vec[p] <- sigma2 / n_rep_post
    } else {
      sigma2_vec[p] <- sigma2 / n_rep_post
    }
  }

  R <- diag(sigma2_vec) + sigma2 / n_rep_pre *
    matrix(1, p, p)

  diag(R) <- sigma2_vec + sigma2 / n_rep_pre

  R_correct <- matrix(0, p, p)
  for (j in seq_len(p)) {
    for (l in seq_len(p)) {
      var_ej <- sigma2_vec[j]
      var_el <- sigma2_vec[l]
      var_e0 <- sigma2 / n_rep_pre
      R_correct[j, l] <- if (j == l) {
        var_ej + var_e0
      } else {
        var_e0
      }
    }
  }

  Ri <- solve(R_correct)

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
