library(tinytest)

# marginal_reduction is positive for favorable params
dr <- marginal_reduction(1, 2, 0, 0.0025, 0.9745, 1.0)
expect_true(dr > 0)


# marginal_reduction eventually decreases
dr3 <- marginal_reduction(3, 2, 0, 1.0, 1.0, 1.0)
dr4 <- marginal_reduction(4, 2, 0, 1.0, 1.0, 1.0)
dr5 <- marginal_reduction(5, 2, 0, 1.0, 1.0, 1.0)
expect_true(dr3 > dr4)
expect_true(dr4 > dr5)


# optimal_allocation returns valid partition
opt <- optimal_allocation(5, 1.0, 1.0, 1.0)
expect_equal(opt$r + opt$k + opt$f, 5)
expect_true(opt$k >= 1)
expect_true(opt$r >= 0)
expect_true(opt$f >= 0)


# optimal_allocation finds minimum variance
opt <- optimal_allocation(4, 1.0, 1.0, 1.0)
v_alt <- var_gamma_matrix(0, 4, 0, 1.0, 1.0, 1.0)
expect_true(opt$variance <= v_alt)


# common_close_beneficial returns logical
result <- common_close_beneficial(2, 2, 1.0, 1.0, 1.0)
expect_true(is.logical(result))

