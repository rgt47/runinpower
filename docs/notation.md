# Notation Reference

Notational conventions used in the manuscript and derived from
Frost, Kenward, and Fox (2008, *Statistics in Medicine*,
27:3717--3731).

## Trial structure

| Symbol | Meaning |
|--------|---------|
| $N$ | Total number of participants |
| $m$ | Number of participants per group |
| $J_0$ | Number of run-in (pre-randomization) observations |
| $J_1$ | Number of post-randomization (treatment period) observations |
| $J_2$ | Number of common close (post-treatment) observations |
| $J$ | Total change scores: $J = J_0 + J_1 + J_2$ |
| $t$ | Equally spaced observation interval |
| $t_j$ | Time of observation $j$: $t_j = jt$ |

## Indices and indicators

| Symbol | Meaning |
|--------|---------|
| $i$ | Participant index |
| $j$ | Time point index (negative for run-in, 0 for baseline, positive for post-randomization) |
| $g_i$ | Treatment group indicator: 1 = treatment, 0 = placebo |
| $h_j$ | Phase indicator: 0 if $j \le 0$ (run-in), 1 if $j > 0$ (post-randomization) |

## Outcome variables

| Symbol | Meaning |
|--------|---------|
| $Y_{ij}$ | Outcome for participant $i$ at time $j$ |
| $c_{ij}$ | Change score: $c_{ij} = Y_{ij} - Y_{i0}$ |
| $C_{ij}$ | Change between consecutive visits (Frost notation): $C_{ij} = Y_{ij} - Y_{i(j-1)}$ |

## Fixed effects

| Symbol | Meaning |
|--------|---------|
| $\alpha$ | Overall intercept (eliminated by change scores) |
| $\delta$ | Rate of change during run-in (common to all participants) |
| $\beta$ | Rate of change during post-randomization for placebo |
| $\gamma$ | Treatment effect on rate of change (primary estimand) |
| $\hat{\gamma}$ | GLS estimator of $\gamma$ |
| $\Delta$ | Minimum clinically important treatment effect (for power calculations) |

## Random effects

| Symbol | Meaning |
|--------|---------|
| $a_i$ | Random intercept: $a_i \sim N(0, \sigma_a^2)$ |
| $b_i$ | Random slope: $b_i \sim N(0, \sigma_b^2)$ |
| $e_{ij}$ | Residual error: $e_{ij} \sim N(0, \sigma^2)$, independent |

## Variance parameters

| Symbol | Meaning |
|--------|---------|
| $\sigma^2$ | Residual (measurement error) variance |
| $\sigma_b^2$ | Between-subject random slope variance |
| $\sigma_a^2$ | Random intercept variance (eliminated by change scores) |
| $\sigma_u^2$ | Within-subject between-visit effects (Frost Table II; not used in the manuscript's simplified model) |
| $a$ | Shorthand for $\sigma^2$ (used in matrix algebra) |
| $b$ | Shorthand for $\sigma_b^2$ (used in matrix algebra) |

## Matrix notation

| Symbol | Meaning |
|--------|---------|
| $X$ | Fixed effects design matrix (columns for $\delta$, $\beta$, $\gamma$; or $\beta$, $\gamma$ when $J_0=0$) |
| $X_i^{(\text{trt})}$, $X_i^{(\text{plc})}$ | Per-subject design matrices for treatment and placebo |
| $Z$ | Random effects design vector: $Z = t \cdot \mathbf{d}$ where $\mathbf{d}$ is the vector of time distances from baseline |
| $G$ | Random effects covariance: $G = \sigma_b^2$ (scalar) |
| $R$ | Residual covariance of change scores: $R = \sigma^2(\mathbf{I}_J + \mathbf{1}_J\mathbf{1}_J')$ |
| $\Sigma$ | Marginal covariance: $\Sigma = R + ZGZ'$ |
| $H$ | Woodbury intermediate: $H = (G^{-1} + Z'R^{-1}Z)^{-1}$ |
| $W$ | Woodbury correction: sum of per-subject outer products $\sum_g X_i^{(g)'} R^{-1}Z \cdot H \cdot Z'R^{-1}X_i^{(g)}$ |
| $\mathbf{I}_J$ | $J \times J$ identity matrix |
| $\mathbf{1}_J$ | $J$-vector of ones |

## Summation shorthands

| Symbol | Meaning |
|--------|---------|
| $S_{J_0}$ | $\sum_{s=1}^{J_0} s^2 = J_0(J_0+1)(2J_0+1)/6$ |
| $S_{J_1}$ | $\sum_{l=1}^{J_1} l^2 = J_1(J_1+1)(2J_1+1)/6$ |
| $D_{J_0}$ | $\sum_{s=1}^{J_0} s = J_0(J_0+1)/2$ |
| $D_{J_1}$ | $\sum_{l=1}^{J_1} l = J_1(J_1+1)/2$ |
| $\eta_J$ | Scaling factor: $\eta_J = 1 + J t^2 b / [a(J+1)]$ |

## Key identities

| Identity | Equation |
|----------|----------|
| Woodbury | $\Sigma^{-1} = R^{-1} - R^{-1}Z(G^{-1} + Z'R^{-1}Z)^{-1}Z'R^{-1}$ |
| Sherman-Morrison (for $R^{-1}$) | $R^{-1} = \frac{1}{\sigma^2}(\mathbf{I}_J - \frac{1}{J+1}\mathbf{1}_J\mathbf{1}_J')$ |
| GLS estimator | $\hat{\boldsymbol{\beta}} = (X'\Sigma^{-1}X)^{-1}X'\Sigma^{-1}\mathbf{c}$ |
| GLS variance | $\text{Var}(\hat{\boldsymbol{\beta}}) = (X'\Sigma^{-1}X)^{-1}$ |
| Information decomposition | $X'\Sigma^{-1}X = X'R^{-1}X - W$ |

## Key results

| Design | $\text{Var}(\hat{\gamma})$ per participant pair |
|--------|-------|
| Standard ($J_0=0$, $J_1=2$, $J_2=0$) | $(\sigma^2 + 2\sigma_b^2 t^2)/t^2$ |
| Two run-in ($J_0=2$, $J_1=1$, $J_2=0$) | $8\sigma^2(\sigma^2 + 5\sigma_b^2 t^2) / [3t^2(\sigma^2 + 2\sigma_b^2 t^2)]$ |
| General ($J_0$, $J_1$, $J_2$) | Computed numerically via Woodbury identity |

## Reparameterizations in Mathematica files

The `.m` files use two derived quantities:

| Symbol | Definition | Context |
|--------|------------|---------|
| $c$ | $\sigma^2 - \sigma_b^2 t^2$ (i.e., $a - bt^2$) | Appears in simplified $\Sigma^{-1}$ |
| $d$ | $2\sigma^2 + \sigma_b^2 t^2$ (i.e., $2a + bt^2$) | Diagonal of marginalized covariance |

These satisfy $a = (d - c)/3$ and $b = (d + 2c)/(3t^2)$.

## Power and sample size

| Symbol | Meaning |
|--------|---------|
| $\alpha$ | Significance level (two-sided; context distinguishes from the intercept) |
| $1 - \beta_{\text{pow}}$ | Statistical power (subscript distinguishes from placebo slope $\beta$) |
| $z_q$ | Standard normal $q$-th quantile |
| $\text{RE}(J_0, J_1, J_2)$ | Relative efficiency vs. standard design: $\text{Var}(\hat{\gamma})\big|_{0, J_1, 0} / \text{Var}(\hat{\gamma})\big|_{J_0, J_1, J_2}$ |

Sample size per group:
$m = (z_{\alpha/2} + z_{1-\beta_{\text{pow}}})^2 \cdot \text{Var}_1(\hat{\gamma}) / \Delta^2$

## MIRIAD parameter values (Frost Table II)

| Parameter | Value | Description |
|-----------|-------|-------------|
| $\beta$ | $-2.2261$ | Mean annual atrophy rate (% per year) |
| $\sigma_b^2$ | $0.9745$ | Between-subject variance in atrophy rates |
| $\sigma_u^2$ | $0.3338$ | Within-subject between-visit effects |
| $\sigma^2$ | $0.0025$ | Additional measurement error |

The manuscript uses $\sigma_b^2 = 0.9745$ and $\sigma^2 = 0.0025$
with $t = 1$ year, giving $\sigma_b/\sigma \approx 19.7$, which
is well above the critical ratio of ~4 where run-in designs
become advantageous.

## Notation change: observation counts

An earlier draft used single lower-case letters for the three
observation counts. The current manuscript uses indexed capitals
to avoid a notational collision with Frost's $r$ (Pearson
correlation) and to make the three counts read as a coherent
family. The mapping is:

| Earlier symbol | Current symbol | Meaning |
|----------------|----------------|---------|
| $r$ | $J_0$ | Run-in observations |
| $k$ | $J_1$ | Post-randomization observations |
| $f$ | $J_2$ | Common-close observations |
| $p$ | $J$ | Total ($J = J_0 + J_1 + J_2$) |
| $S_r$, $S_k$ | $S_{J_0}$, $S_{J_1}$ | Sum-of-squares shorthands |
| $D_r$, $D_k$ | $D_{J_0}$, $D_{J_1}$ | Triangular-number shorthands |
| $\eta_p$ | $\eta_J$ | Scaling factor |
| $\mathrm{RE}(r,k,f)$ | $\mathrm{RE}(J_0, J_1, J_2)$ | Relative efficiency |

R function arguments in the implementation (e.g.
`var_gamma_matrix(r, k, f)`) retain the lower-case names for
backward compatibility; the rename applies only to mathematical
and narrative notation in the manuscripts.
