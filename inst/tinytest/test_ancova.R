library(tinytest)

# ANCOVA variance is between GLS and averaged
a <- 1.0
b <- 1.0
tt <- 1.0
v_gls <- var_gamma_matrix(2, 2, 0, a, b, tt)
v_anc <- var_gamma_ancova(2, 2, 0, a, b, tt)
v_avg <- var_gamma_avg(2, 2, 0, a, b, tt)
expect_true(v_gls <= v_anc + 1e-10)
expect_true(v_anc <= v_avg + 1e-10)


# ANCOVA equals averaged when r = 0
v_anc <- var_gamma_ancova(0, 2, 0, 1.0, 0.5, 1.0)
v_avg <- var_gamma_avg(0, 2, 0, 1.0, 0.5, 1.0)
expect_equal(v_anc, v_avg)


# re_ancova_vs_gls is <= 1
re <- re_ancova_vs_gls(2, 2, 0, 1.0, 1.0, 1.0)
expect_true(re <= 1 + 1e-10)


# ANCOVA with run-in improves over no run-in
v_no_ri <- var_gamma_ancova(0, 2, 0, 1.0, 1.0, 1.0)
v_ri <- var_gamma_ancova(2, 2, 0, 1.0, 1.0, 1.0)
expect_true(v_ri < v_no_ri)


# ANCOVA is finite for MIRIAD parameters
v <- var_gamma_ancova(2, 2, 0, 0.0025, 0.9745, 1.0)
expect_true(is.finite(v))
expect_true(v > 0)


# ANCOVA works with common close
v <- var_gamma_ancova(2, 2, 2, 1.0, 1.0, 1.0)
expect_true(is.finite(v))
expect_true(v > 0)

