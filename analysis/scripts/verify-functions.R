## Independent numerical verification of runinpower functions.
## Each check builds the quantity from first principles and compares
## to the package. PASS if relative error < tol.
suppressMessages(library(runinpower))
set.seed(20260623)
tol <- 1e-8
pass <- function(name, a, b, t = tol) {
  re <- abs(a - b) / max(abs(b), 1e-12)
  cat(sprintf('%-52s %-9s  pkg=%.6g  indep=%.6g  relerr=%.1e\n',
              name, if (re < t) 'PASS' else '**FAIL**', a, b, re))
}

## Independent per-pair GLS variance: full Sigma inverse, no Woodbury.
gls_bruteforce <- function(Xt, Xp, Sigma, nc) {
  Si <- solve(Sigma)
  info <- t(Xt) %*% Si %*% Xt + t(Xp) %*% Si %*% Xp
  solve(info)[nc, nc]
}

cat('==== 1. var_gamma_matrix vs brute-force full-inverse GLS ====\n')
for (cfg in list(c(0,2,0), c(2,1,0), c(2,2,0), c(3,2,1), c(4,3,2), c(1,2,0))) {
  r <- cfg[1]; k <- cfg[2]; f <- cfg[3]
  a <- 1.3; b <- 0.7; tt <- 0.9
  des <- build_design(r, k, f, tt)
  Sigma <- build_R_matrix(r+k+f, a) + des$Z %*% (b * t(des$Z))
  indep <- gls_bruteforce(des$X_trt, des$X_plc, Sigma, des$ncol_X)
  pkg <- var_gamma_matrix(r, k, f, a, b, tt)
  pass(sprintf('  matrix (r=%d,k=%d,f=%d)', r, k, f), pkg, indep)
}

cat('==== 2. Closed forms vs matrix engine ====\n')
a <- 1.0; b <- 0.5; tt <- 1.2
pass('  Frost closed = matrix(0,2,0)',
     var_gamma_frost(a, b, tt), var_gamma_matrix(0, 2, 0, a, b, tt))
pass('  r2 closed = matrix(2,1,0)',
     var_gamma_r2(a, b, tt), var_gamma_matrix(2, 1, 0, a, b, tt))
## independent algebra of the two closed forms
pass('  Frost closed = (a+2bt^2)/t^2',
     var_gamma_frost(a, b, tt), (a + 2*b*tt^2)/tt^2)
pass('  r2 closed = 8a(a+5bt^2)/(3t^2(a+2bt^2))',
     var_gamma_r2(a, b, tt),
     8*a*(a+5*b*tt^2)/(3*tt^2*(a+2*b*tt^2)))

cat('==== 3. var_gamma_ar1: independent AR(1) raw-process build ====\n')
ar1_change_cov <- function(r, k, f, sigma2, rho) {
  ## raw AR(1) over visit indices incl baseline 0; difference out 0
  idx <- c(if (r>0) -r:-1, if (k>0) 1:k, if (f>0) (k+1):(k+f))
  all_t <- c(0, idx)
  Sraw <- outer(all_t, all_t, function(u, v) sigma2 * rho^abs(u - v))
  L <- cbind(-1, diag(length(idx)))          # c_j = e_idx_j - e_0
  L %*% Sraw %*% t(L)
}
for (cfg in list(c(2,2,0), c(4,2,0), c(2,2,2), c(0,2,0))) {
  r <- cfg[1]; k <- cfg[2]; f <- cfg[3]
  a <- 1.1; b <- 0.6; tt <- 1.0; rho <- 0.5
  Rc <- ar1_change_cov(r, k, f, a, rho)
  des <- build_design(r, k, f, tt)
  Sigma <- Rc + des$Z %*% (b * t(des$Z))
  indep <- gls_bruteforce(des$X_trt, des$X_plc, Sigma, des$ncol_X)
  pkg <- var_gamma_ar1(r, k, f, a, b, tt, rho)
  pass(sprintf('  ar1 (r=%d,k=%d,f=%d,rho=.5)', r, k, f), pkg, indep)
}
## AR(1) at rho->0 should approach compound-symmetry matrix result
pass('  ar1(rho~0) ~ CS matrix',
     var_gamma_ar1(2, 2, 0, 1, .5, 1, 1e-6),
     var_gamma_matrix(2, 2, 0, 1, .5, 1), t = 1e-3)

cat('==== 4. var_gamma_avg vs Monte Carlo (raw simulation) ====\n')
mc_avg <- function(r, k, f, a, b, tt, nrep = 80000, m = 400) {
  des <- build_design(r, k, f, tt); p <- r+k+f
  Sigma <- build_R_matrix(p, a) + des$Z %*% (b * t(des$Z))
  Lc <- chol(Sigma)
  w <- if (r == 0) rep(1/(k+f), p) else
    c(rep(-1/r, r), rep(1/(k+f), k+f))
  ## per-subject averaged change score d = w' c; var across subjects
  ds <- w %*% (t(Lc) %*% matrix(rnorm(p*nrep), p))
  ## Var(gamma_avg) = Var(dbar_trt - dbar_plc) for m per grp = 2 Var(d)/m,
  ## but package returns per-pair (m=1) scale: 2 Var(d). Compare variances.
  2 * var(as.numeric(ds))
}
for (cfg in list(c(2,2,0), c(2,2,2), c(0,2,0))) {
  r <- cfg[1]; k <- cfg[2]; f <- cfg[3]
  a <- 1.0; b <- 1.0; tt <- 1.0
  pass(sprintf('  avg (r=%d,k=%d,f=%d) MC', r, k, f),
       var_gamma_avg(r, k, f, a, b, tt),
       mc_avg(r, k, f, a, b, tt), t = 3e-2)
}

cat('==== 5. var_gamma_ancova vs Monte Carlo ANCOVA regression ====\n')
mc_ancova <- function(r, k, f, a, b, tt, nrep = 4000, m = 4000) {
  des <- build_design(r, k, f, tt); p <- r+k+f
  Sigma <- build_R_matrix(p, a) + des$Z %*% (b * t(des$Z))
  Lc <- chol(Sigma)
  est <- numeric(nrep)
  for (s in seq_len(nrep)) {
    Ct <- t(Lc) %*% matrix(rnorm(p*m), p)   # change scores, trt
    Cp <- t(Lc) %*% matrix(rnorm(p*m), p)   # placebo
    post <- function(C) colMeans(C[(r+1):p, , drop=FALSE])
    pre  <- function(C) if (r==0) rep(0, ncol(C)) else
      colMeans(C[1:r, , drop=FALSE])
    y <- c(post(Ct), post(Cp)); g <- c(rep(1,m), rep(0,m))
    if (r == 0) {
      est[s] <- coef(lm(y ~ g))['g']
    } else {
      xpre <- c(pre(Ct), pre(Cp))
      est[s] <- coef(lm(y ~ g + xpre))['g']
    }
  }
  var(est) * m   # scale to per-pair (m per group) -> *m gives per-pair
}
for (cfg in list(c(2,2,0), c(0,2,0))) {
  r <- cfg[1]; k <- cfg[2]; f <- cfg[3]
  a <- 1.0; b <- 1.0; tt <- 1.0
  pass(sprintf('  ancova (r=%d,k=%d,f=%d) MC', r, k, f),
       var_gamma_ancova(r, k, f, a, b, tt),
       mc_ancova(r, k, f, a, b, tt), t = 5e-2)
}

cat('==== 6. Efficiency ordering GLS <= ANCOVA <= averaged ====\n')
ord_ok <- TRUE
for (r in c(0,1,2,4)) for (k in c(1,2,3)) for (f in c(0,1)) {
  if (r+k+f < 2) next
  for (ratio in c(0.5,1,2,5,20)) {
    b <- ratio^2
    vg <- var_gamma_matrix(r,k,f,1,b,1)
    va <- var_gamma_ancova(r,k,f,1,b,1)
    vv <- var_gamma_avg(r,k,f,1,b,1)
    if (vg > va*(1+1e-9) || va > vv*(1+1e-9)) {
      ord_ok <- FALSE
      cat(sprintf('   ORDER VIOLATION r=%d k=%d f=%d ratio=%g: GLS=%.4f ANC=%.4f AVG=%.4f\n',
                  r,k,f,ratio,vg,va,vv))
    }
  }
}
cat(sprintf('  ordering holds across grid: %s\n', if (ord_ok) 'PASS' else '**FAIL**'))

cat('==== 7. var_gamma_replicated sanity ====\n')
pass('  replicated(1,1) = matrix',
     var_gamma_replicated(0,2,0,1,1,1,n_rep_post=1,n_rep_pre=1),
     var_gamma_matrix(0,2,0,1,1,1))
cat(sprintf('  monotone in n_rep_post (1>=2>=5)? %s\n',
  {v1<-var_gamma_replicated(0,2,0,1,0.2,1,n_rep_post=1)
   v2<-var_gamma_replicated(0,2,0,1,0.2,1,n_rep_post=2)
   v5<-var_gamma_replicated(0,2,0,1,0.2,1,n_rep_post=5)
   if (v1>=v2 && v2>=v5) 'PASS' else '**FAIL**'}))

cat('==== 8. Paper 2 chunk ar1_var_mean: suspected off-by-c ====\n')
## chunk-as-written:
ar1_var_mean_chunk <- function(n, sigma2, rho) {
  if (n == 1) return(sigma2)
  c <- rho
  S_n <- (n*(1-c) - (1-c^n))/(1-c)^2
  sigma2/n * (1 + 2*rho/(n*c) * S_n)
}
## independent truth: Var(mean of n consecutive AR(1)) =
## sigma2/n^2 * sum_{i,j} rho^|i-j|
ar1_var_mean_true <- function(n, sigma2, rho) {
  M <- outer(1:n, 1:n, function(i,j) rho^abs(i-j))
  sigma2 / n^2 * sum(M)
}
for (n in c(2,3,5)) for (rho in c(0.3,0.5,0.7)) {
  pass(sprintf('  ar1_var_mean n=%d rho=%.1f', n, rho),
       ar1_var_mean_chunk(n, 1, rho),
       ar1_var_mean_true(n, 1, rho), t = 1e-6)
}
cat('\nDONE\n')
