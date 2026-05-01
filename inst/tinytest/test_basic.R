library(tinytest)

# build_R_matrix produces correct structure
R <- build_R_matrix(3, 1.0)
expect_equal(dim(R), c(3, 3))
expect_equal(diag(R), rep(2, 3))
expect_equal(R[1, 2], 1.0)
expect_equal(R[1, 3], 1.0)


# build_design returns correct dimensions
des <- build_design(2, 1, 0, 1.0)
expect_equal(des$p, 3)
expect_equal(des$ncol_X, 3)
expect_equal(nrow(des$X_trt), 3)
expect_equal(nrow(des$X_plc), 3)
expect_equal(nrow(des$Z), 3)


# build_design omits delta when r = 0
des <- build_design(0, 2, 0, 1.0)
expect_equal(des$ncol_X, 2)
expect_equal(colnames(des$X_trt), c("beta", "gamma"))


# var_gamma_matrix matches Frost closed form
a <- 1.0
b <- 0.5
tt <- 1.0
v_matrix <- var_gamma_matrix(0, 2, 0, a, b, tt)
v_closed <- var_gamma_frost(a, b, tt)
expect_equal(v_matrix, v_closed)


# var_gamma_matrix matches r=2 closed form
a <- 1.0
b <- 0.5
tt <- 1.0
v_matrix <- var_gamma_matrix(2, 1, 0, a, b, tt)
v_closed <- var_gamma_r2(a, b, tt)
expect_equal(v_matrix, v_closed)


# run-in reduces variance when sigma_b/sigma is large
a <- 0.0025
b <- 0.9745
tt <- 1.0
v0 <- var_gamma_matrix(0, 2, 0, a, b, tt)
v2 <- var_gamma_matrix(2, 2, 0, a, b, tt)
expect_true(v2 < v0)


# relative_efficiency > 1 for favorable parameters
re <- relative_efficiency(2, 2, 0, 0.0025, 0.9745, 1.0)
expect_true(re > 1)


# sample_size returns positive integer
n <- sample_size(0, 2, 0, 1.0, 0.5, 1.0, 0.5)
expect_true(is.numeric(n))
expect_true(n > 0)
expect_equal(n, as.integer(n))


# power_calc is monotone in m
pow10 <- power_calc(0, 2, 0, 1.0, 0.5, 1.0, 0.5, 10)
pow50 <- power_calc(0, 2, 0, 1.0, 0.5, 1.0, 0.5, 50)
expect_true(pow50 > pow10)

