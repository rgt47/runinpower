# Independent Numerical Verification of `runinpower`
*2026-06-23 19:20 PDT*

Each package quantity was rebuilt from first principles and compared
to the package output. Script: reproducible with `set.seed(20260623)`.

## Verified correct (PASS)

| Function | Independent check | Result |
|---|---|---|
| `var_gamma_matrix` | full-Sigma-inverse GLS (no Woodbury), 6 configs incl. common close | exact (relerr < 1e-15) |
| `var_gamma_frost` | matrix engine and `(sigma^2+2 sigma_b^2 t^2)/t^2` | exact |
| `var_gamma_r2` | matrix engine and `8a(a+5bt^2)/(3t^2(a+2bt^2))` | exact (re-derived by hand earlier) |
| `var_gamma_ar1` | raw AR(1) process built over visit indices, differenced to change scores, full GLS; 4 configs | exact; `rho->0` recovers the CS result |
| `build_R_matrix_ar1` | equals `L Sigma_raw L'` differencing construction | exact |
| `var_gamma_avg` | Monte Carlo over simulated change scores | match within MC error (relerr < 4e-3) |
| `var_gamma_ancova` | Monte Carlo ANCOVA regression `lm(post ~ g + pre)` | match within MC error (relerr < 2e-2) |
| GLS <= ANCOVA <= averaged | full grid (r in 0,1,2,4; k in 1,2,3; f in 0,1; sigma_b/sigma in 0.5..20) | ordering holds everywhere |
| `var_gamma_replicated` | equals `var_gamma_matrix` at `n_rep=1`; monotone decreasing in `n_rep` | PASS |

The Paper 1-4 numerical tables and figures all derive from
`var_gamma_matrix` / `var_gamma_ar1` (and the closed forms), so the
substantive results of all four papers are numerically sound.

## Bug found and fixed (Paper 2, moment-based section only)

The illustrative `moment_based` chunk in Paper 2 (Section 4) is the one
place with errors. It does not feed any table or figure (it prints one
`cat()` comparison), so no headline result is affected, but it was wrong
as written:

1. **`ar1_var_mean` off by a factor.** As written,
   `sigma2/n * (1 + 2*rho/(n*c) * S_n)` overstates the variance of an
   AR(1) phase mean (n=2, rho=0.5 gave 1.0; truth `(1+rho)/2 = 0.75`).
   The independent truth is `sigma2/n^2 * sum_{i,j} rho^|i-j|`. The
   error is the spurious `/c`. **Fixed** to `2*rho/n * S_n`, which now
   matches the truth exactly for all tested `(n, rho)`.

2. **`ar1_cov_phase_means` called with the wrong `K`.** The closed-form
   cross-covariance is correct only when `K` is the time-index lag from
   the last run-in visit (-1) to the first treatment visit (+1), which
   is **2** in this model (baseline sits at 0 between them). The chunk
   passed `K=1`, off by one, producing a value too small by `1/rho`.
   **Fixed** the call to `K=2`. (The displayed cross-covariance formula
   was separately aligned to the code in the earlier consistency pass.)

These fixes make the moment-based illustration internally correct. The
section as a whole remains a draft (its general `(J_0,J_1,J_2)`
derivation is still "to be derived"); a full re-derivation is left to
the author.

## Not independently re-derived

The sandwich-variance functions in Paper 4's misspecification chunks
(`var_misspec_gls`, `var_misspec_r`) implement the standard
`bread %*% meat %*% bread` form and were read but not checked against an
external sandwich computation. The `optimal_allocation` and
`common_close_beneficial` routines are thin brute-force wrappers over
the verified `var_gamma_matrix` and require no separate check.

---
*Rendered on 2026-06-23 at 19:20 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/analysis/report/verification-findings-2026-06-23.md*
