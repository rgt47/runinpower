library(tinytest)

# averaged variance is >= GLS variance
v_gls <- var_gamma_matrix(2, 2, 0, 1.0, 1.0, 1.0)
v_avg <- var_gamma_avg(2, 2, 0, 1.0, 1.0, 1.0)
expect_true(v_avg >= v_gls - 1e-10)


# re_avg_vs_gls is <= 1
re <- re_avg_vs_gls(2, 2, 0, 1.0, 1.0, 1.0)
expect_true(re <= 1 + 1e-10)


# averaged variance is finite and positive
v <- var_gamma_avg(2, 2, 0, 0.0025, 0.9745, 1.0)
expect_true(is.finite(v))
expect_true(v > 0)


# averaged variance works without run-in
v <- var_gamma_avg(0, 2, 0, 1.0, 0.5, 1.0)
expect_true(is.finite(v))
expect_true(v > 0)


# averaged variance works with common close
v <- var_gamma_avg(2, 2, 2, 1.0, 1.0, 1.0)
expect_true(is.finite(v))
expect_true(v > 0)

