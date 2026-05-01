library(tinytest)

# no replication matches standard variance
v_std <- var_gamma_matrix(2, 2, 0, 1.0, 1.0, 1.0)
v_rep <- var_gamma_replicated(2, 2, 0, 1.0, 1.0, 1.0,
                                 n_rep_post = 1,
                                 n_rep_pre = 1)
expect_equal(v_rep, v_std, tolerance = 1e-8)


# endpoint replication reduces variance
v1 <- var_gamma_replicated(0, 2, 0, 1.0, 0.5, 1.0,
                              n_rep_post = 1)
v3 <- var_gamma_replicated(0, 2, 0, 1.0, 0.5, 1.0,
                              n_rep_post = 3)
expect_true(v3 < v1)


# baseline replication reduces variance
v1 <- var_gamma_replicated(0, 2, 0, 1.0, 0.5, 1.0,
                              n_rep_pre = 1)
v3 <- var_gamma_replicated(0, 2, 0, 1.0, 0.5, 1.0,
                              n_rep_pre = 3)
expect_true(v3 < v1)


# replication benefit is larger when sigma2 dominates
re_high_noise <- var_gamma_replicated(
    0, 2, 0, 5.0, 0.1, 1.0, n_rep_post = 3) /
    var_gamma_replicated(
      0, 2, 0, 5.0, 0.1, 1.0, n_rep_post = 1)

re_low_noise <- var_gamma_replicated(
    0, 2, 0, 0.1, 5.0, 1.0, n_rep_post = 3) /
    var_gamma_replicated(
      0, 2, 0, 0.1, 5.0, 1.0, n_rep_post = 1)

expect_true(re_high_noise < re_low_noise)


# replicated variance is finite and positive
v <- var_gamma_replicated(2, 2, 0, 0.0025, 0.9745,
                             1.0, n_rep_post = 3,
                             n_rep_pre = 2)
expect_true(is.finite(v))
expect_true(v > 0)

