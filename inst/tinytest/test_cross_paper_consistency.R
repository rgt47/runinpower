library(tinytest)

# Cross-paper consistency checks.
#
# The four manuscripts (analysis/report/{01,02,03,04}-*/) share a
# single computational engine built around var_gamma_matrix()
# (Paper 1). This file asserts that the engines used by Papers 2-4
# agree with Paper 1's Woodbury routine at parameter points where
# the underlying math implies exact agreement, guarding against a
# future refactor silently desynchronizing one paper's numbers from
# the shared core. Added per pub_review whitepaper item (c)10,
# 2026-08-16 remediation.

a  <- 1.0
b  <- 0.5
tt <- 1.0

# Paper 2 (AR(1)) vs Paper 1 (Woodbury/compound-symmetry):
# at rho = 0, build_R_matrix_ar1() reduces algebraically to the
# compound-symmetry residual covariance sigma2 * (I + 11'), the
# same R used by var_gamma_matrix(); the two routines must therefore
# agree exactly at rho = 0 for every (r, k, f) configuration.
for (cfg in list(c(0, 2, 0), c(2, 1, 0), c(2, 2, 2), c(4, 3, 1))) {
  r <- cfg[1]; k <- cfg[2]; f <- cfg[3]
  v_matrix <- var_gamma_matrix(r, k, f, a, b, tt)
  v_ar1_rho0 <- var_gamma_ar1(r, k, f, a, b, tt, rho = 0)
  expect_equal(
    v_ar1_rho0, v_matrix, tolerance = 1e-10,
    info = sprintf(
      "AR(1) engine at rho=0 must match Woodbury engine at (%d,%d,%d)",
      r, k, f))
}

# Paper 3 (allocation) vs Paper 1 (Woodbury): marginal_reduction()
# is defined as a direct difference of two var_gamma_matrix() calls,
# so it must reproduce that difference exactly for any r >= 1.
r <- 2; k <- 2; f <- 0
expect_equal(
  marginal_reduction(r, k, f, a, b, tt),
  var_gamma_matrix(r - 1, k, f, a, b, tt) -
    var_gamma_matrix(r, k, f, a, b, tt),
  tolerance = 1e-12,
  info = "allocation engine must reproduce Woodbury variance differences")

# Paper 4 (ANCOVA/GLS relative efficiency) vs Paper 1 (Woodbury):
# re_ancova_vs_gls() is defined as var_gamma_matrix() /
# var_gamma_ancova(), so it must reproduce that ratio exactly, and
# ANCOVA can never be more efficient than GLS (RE <= 1) at a shared
# parameter point, per the ordering proposition in Paper 4.
r <- 2; k <- 2; f <- 0
re <- re_ancova_vs_gls(r, k, f, a, b, tt)
expect_equal(
  re,
  var_gamma_matrix(r, k, f, a, b, tt) /
    var_gamma_ancova(r, k, f, a, b, tt),
  tolerance = 1e-12,
  info = "Paper 4 RE must reproduce the Woodbury/ANCOVA variance ratio")
expect_true(
  re <= 1 + 1e-8,
  info = "GLS must be at least as efficient as ANCOVA (RE <= 1)")
