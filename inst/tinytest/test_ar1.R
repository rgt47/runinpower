library(tinytest)

# AR(1) covariance reduces to CS when rho = 0
d_all <- c(-2, -1, 1)
R_ar1 <- build_R_matrix_ar1(3, 1.0, 0, d_all)
R_cs <- build_R_matrix(3, 1.0)
expect_equal(R_ar1, R_cs)


# AR(1) covariance diagonal is correct
d_all <- c(1, 2)
R <- build_R_matrix_ar1(2, 1.0, 0.5, d_all)
expect_equal(R[1, 1], 1 * (1 - 0.5 - 0.5 + 1))
expect_equal(R[2, 2], 1 * (1 - 0.5^2 - 0.5^2 + 1))


# AR(1) variance matches CS at rho = 0
a <- 1.0
b <- 0.5
tt <- 1.0
v_cs <- var_gamma_matrix(2, 2, 0, a, b, tt)
v_ar1 <- var_gamma_ar1(2, 2, 0, a, b, tt, rho = 0)
expect_equal(v_ar1, v_cs, tolerance = 1e-10)


# AR(1) variance changes with rho
a <- 1.0
b <- 0.5
tt <- 1.0
v_low <- var_gamma_ar1(2, 2, 0, a, b, tt, rho = 0.1)
v_high <- var_gamma_ar1(2, 2, 0, a, b, tt, rho = 0.9)
expect_false(isTRUE(all.equal(v_low, v_high)))

