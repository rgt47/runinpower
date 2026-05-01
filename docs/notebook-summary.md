# Summary of Notebook Derivations

This document indexes the distinct mathematical developments
found across six scanned notebooks (2019--2024) and their
LaTeX transcriptions. Repeated derivations are noted once
at their most complete appearance.

## 1. Foundational framework (Aug--Oct 2019)

**Transcribed in:** `notebook-foundations.tex`

### 1.1 General linear mixed model

The raw-outcome model with random intercept and slope:

$$Y_{ij} = \alpha + \alpha_{0i} + (\beta_1 + \beta_{1i}
  + \beta_x g_i) t_j + e_{ij}$$

Matrix form $Y_i = X_i\beta + Z_i\gamma_i + \epsilon_i$,
marginal covariance $\Sigma = R + ZDZ'$.

*Source: 8/31/19, 9/6/19.*

### 1.2 Variance with both random intercept and slope

General formula including $\sigma_\alpha^2$, $\sigma_\beta^2$,
and their correlation $\rho$:

$$\text{Var}(\hat\beta_{TX}) = \frac{2}{m}\left[
  \frac{\sigma^2/K + \sigma_\beta^2 + \sigma_\alpha^2
    + \frac{\sigma_\alpha^2\sigma_\beta^2}{\sigma^2}(1-\rho^2)}
  {1 + K\sigma_\alpha^2/\sigma^2}\right]$$

Reduces to $\frac{2}{m}[\sigma^2/K + \sigma_\beta^2]$
when $\sigma_\alpha^2 = 0$.

*Source: 9/6/19.*

### 1.3 Change-score parameterization

Elimination of intercepts via $z_{ij} = y_{ij} - y_{i0}$,
yielding $z_i = X\beta + Ub_i + \varepsilon_i$ with
$\varepsilon_j \sim N(0, \sigma^2 I)$.

*Source: 9/25/19.*

### 1.4 Seven-step algorithm

Computational recipe for $\text{Var}(\hat\gamma)$:
(1) $V_i^* = R_i + Z_i^*DZ_i^{*'}$,
(2) $V_i^{*-1}$,
(3-4) $X_i'V_i^{-1}X_i$ per group,
(5) sum,
(6) invert,
(7) extract $(2,2)$ element.

*Source: 9/25/19.*

### 1.5 Time normalization

$t_i^* = (t_i - \bar{t})/\text{SD}(t_i)$ with
$\sum t_i^* = 0$, $\sum t_i^{*2} = K-1$.
Original-scale recovery:
$V(\hat\alpha_x) = V(\hat\alpha_x^*)/s_x^2$.

*Source: 9/25/19.*

### 1.6 AR(1) covariance sums

Closed-form for $S_n = \sum_{i=1}^{n-1}(n-i)c^{i-1}$:

$$S_n = \frac{n(1-c) - (1-c^n)}{(1-c)^2}$$

Verified for $n = 2, 3, 4, 5$.

*Source: 10/14/19. Also repeated with additional detail
in scans 4 and 5.*


## 2. Woodbury identity and Frost reproduction (Oct 2019 -- Nov 2021)

**Transcribed in:** `notebook-derivations.tex`

### 2.1 Woodbury inversion of raw-outcome covariance

For $V = \sigma^2 I + \sigma_b^2 \mathbf{t}\mathbf{t}'$
(diagonal $R$):

$$V^{-1} = \frac{1}{a}\left(I
  - \frac{b}{a+b}\mathbf{t}\mathbf{t}'\right)$$

and the information matrix decomposition
$X'V^{-1}X = \frac{1}{a}(X'X - \rho X'\mathbf{t}\mathbf{t}'X)$
where $\rho = b/(a+b)$.

*Source: 10/7/19 (handwritten, scan 1).*

### 2.2 Frost 3.1.1 by hand (conventional design, $r=0$, $k=2$)

Change-score version with $R = a\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$.
Complete Woodbury yielding:

$$\text{Var}(\hat\gamma)\big|_{r=0,k=2}
  = \frac{\sigma^2 + 2\sigma_b^2 t^2}{t^2}$$

Multiple independent derivations exist (10/20/19 "Take 2";
10/21; 8/9/24). The 8/9/24 version also derives the
explicit GLS estimator:

$$\hat\gamma = \frac{1}{2t}
  [(C_{22}+C_{21}) - (C_{12}+C_{11})]$$

*Source: "Take 2" pages (scan 3), 10/21 (scan 2),
8/9/24 (scan 6).*

### 2.3 Frost 3.1.2 by hand (run-in, $r=2$, $k=1$)

Three-parameter model ($\delta, \beta, \gamma$) with
$p = 3$ change scores, $R = a(I_3 + \mathbf{1}_3\mathbf{1}_3')$,
$Z = t(-2, -1, 1)'$. Complete Woodbury:

$$\text{Var}(\hat\gamma)\big|_{r=2,k=1}
  = \frac{8\sigma^2(\sigma^2 + 5\sigma_b^2 t^2)}
  {3t^2(\sigma^2 + 2\sigma_b^2 t^2)}$$

Includes verification of the correlation factor
$r = (b-a)/(b+2a)$ with $1-r^2 = 3a(a+2b)/(b+2a)^2$,
confirmed by checking $(2b+4a)(1-r^2) \cdot 3a/(b+2a) = 6a$.

*Source: 11/11--11/12 (scan 2).*

### 2.4 Reparameterization ($c$, $d$ variables)

$c = a - bt^2$ (or $c = 2a + bt^2$ in some pages),
$d = 2a + bt^2$ (or $d = a + 2bt^2$). These appear in the
Mathematica files and in the 11/16 notes (scan 3).
The manuscript settles on the change-score formulation
without reparameterization.

*Source: frost312.m; 11/16 notes.*

### 2.5 Test statistic and relative sample size

$T(d,r) = \bar{C}_{1q}(d,r) - \bar{C}_{0q}(d,r)$,
relative sample size $N_{T(d,r)}/N_{T(0,0)}$,
sample size formula
$N = (z_{\alpha/2}+z_\beta)^2 V/\Delta^2$.

*Source: 3/10/21 (scan 2).*


## 3. Extensions and generalizations

**Transcribed in:** `notebook-extensions.tex`

### 3.1 General time points ($t_1, t_2$ not equally spaced)

With $\sum t_i = 0$, $\sum t_i^2 = 1$:

$$\text{Var}(\hat\gamma) = \frac{3a}{1 - t_1 t_2} + 2b$$

Woodbury scalar: $H = 3ab/(3a + 2b(1-t_1t_2))$.

*Source: 11/16 (scan 3).*

### 3.2 Conventional design with $K$ visits (raw outcomes)

$$\text{Var}(\hat\gamma) = \frac{2\sigma^2}{m(K-1)}
  (1 + (K-1)\rho)$$

where $\rho = \sigma_b^2/(\sigma^2 + \sigma_b^2)$.
Treatment information: $\frac{K-1}{\sigma^2}
\cdot\frac{1}{1+(K-1)\rho}\begin{pmatrix}1&1\\1&1\end{pmatrix}$.

*Source: scan 3 (TRT/PLA pages); scan 4 (10/4 pages).*

### 3.3 AR(1) cross-covariance between run-in and treatment

$$\text{Cov}(\bar{X}_0, \bar{X}_K)
  = \frac{\sigma^2}{rc}\rho_1 c^{K-1}
  \frac{(1-c^r)(1-c^q)}{(1-c)^2}$$

Uses the Knuth identity (pg 40) for double geometric sums.

*Source: 10/13 (scan 3).*

### 3.4 Three run-in observations ($r=3$, $k=1$)

$p = 4$ change scores. Sherman-Morrison for $R^{-1}$ with
$p+1 = 5$: $R^{-1} = \frac{1}{5a}(5I - \mathbf{1}\mathbf{1}')$.
General information matrix with sums $S_1, \ldots, S_4$.

*Source: 11/12 (scan 5, page 1).*

### 3.5 General $r$ run-in: information matrix with $S_i$ sums

Design matrix with time sums $S_1, S_2, S_3$ (run-in)
and $S_4$ (treatment). Using $\sum S_i = 0$ and
$e = r + K$:

$$X'R^{-1}X = \frac{1}{a(e+1)}\begin{pmatrix}
  S_1^2 e - S_2 - S_3 & S_4^2 & S_4^2 \\
  S_4^2 & S_4^2 e & S_4^2 e \\
  S_4^2 & S_4^2 e & S_4^2 e
\end{pmatrix}$$

This is the precursor to the general formula in the
report's Section 3.3.

*Source: scan 5, page 2.*

### 3.6 Simplest run-in case ($r=1$, $k=1$)

Two change scores (one run-in, one post-randomization) with
the consecutive-change $R$. Comparison of design A
($r=1$, $k=1$) vs B ($r=0$, $k=1$) using the
2-parameter model.

*Source: 11/13/19 and 11/16 (scan 5).*

### 3.7 Simplest Frost case ($r=0$, $k=1$: one baseline, one follow-up)

Single change score $C_1 = Y_1 - Y_0$, scalar $R = 2a$.
$H = ab/(a+bt^2)$.

*Source: scan 5, page 6.*

### 3.8 Stratified variance with unequal visit patterns

Participants have varying $(r_i, f_i)$. Variance computed
by summing within strata weighted by group sizes
$(m_1, m_2, \ldots)$.

*Source: 10/14/19 (scan 3 and notebook-foundations).*


## 4. Explicit correlation structures

### 4.1 Marginal covariance $V = R + ZGZ'$ (computed explicitly)

**Consecutive-change parameterization**
($C_1 = Y_0 - Y_{-2}$, $C_2 = Y_0 - Y_{-1}$,
$C_3 = Y_1 - Y_0$):

$$\text{Corr}(C_1, C_2) = \frac{a + 2t^2 b}{2a + 4t^2 b},
\quad
\text{Corr}(C_1, C_3) = \frac{-a + 2t^2 b}{2a + 4t^2 b},
\quad
\text{Corr}(C_2, C_3) = \frac{-a + t^2 b}{2a + t^2 b}$$

**Change-from-baseline parameterization**
($C_1 = Y_{-2} - Y_0$, $C_2 = Y_{-1} - Y_0$,
$C_3 = Y_1 - Y_0$):

$$\text{Corr}(C_1, C_4) = \frac{-a + bt^2}{2a + bt^2},
\quad
\text{Corr}(C_1, C_3) = \frac{bt^2}{2a + bt^2}$$

*Source: scan 5, pages 7--8.*

### 4.2 Five correlation cases

1. Compound symmetry: $\text{Corr} = \rho$ for all $k \ne k'$
2. AR(1): $\text{Corr} = \rho^{|k'-k|}$
3. Case 1 with variable $r_i$
4. Case 2 with variable $r_i, d_i$
5. General: $\text{Corr} = \sigma_{|k-k'|}/\sigma^2$

*Source: 3/10/21 (scan 2); 5/13/20 and 12/29 (scan 4).*


## 5. GLS estimator (not just variance)

### 5.1 Full GLS estimator for Frost 3.1.1

$$\hat\gamma = \frac{1}{2t}
  [(C_{22} + C_{21}) - (C_{12} + C_{11})]$$

Derived via $(X'\Sigma^{-1}X)^{-1}X'\Sigma^{-1}Y$,
showing the estimator is a contrast of summed
consecutive changes between groups.

*Source: 8/9/24 (scan 6); also frost311.m.*

### 5.2 Full GLS estimator in normalized time

$$\hat\beta_x = \frac{m}{a+b}
  \left(\sum^T t_i y_i - \sum^P t_i y_i\right)$$

Treatment effect is the difference in time-weighted
outcome sums, scaled by $(a+b)^{-1}$.

*Source: 10/6/19 (scan 4).*


## 6. Alternative approaches (not in the manuscript)

### 6.1 SUR (Stanek) formulation

Seemingly Unrelated Regression with change and baseline
as separate equations. Deviation-from-mean
parameterization:
$\begin{pmatrix}\beta_{11}\\\beta_{12}\end{pmatrix}
= \frac{1}{2}\begin{pmatrix}\delta_2+\delta_1\\
\delta_2-\delta_1\end{pmatrix}$.

*Source: 10/27 (scan 4).*

### 6.2 Crossover (N-of-1) design

A/B crossover with sequences A:B and B:A. Model:
$Y_{111} = \mu + \epsilon_{11} + \Pi_1 + \phi_1$.
Full OLS design matrix with 9 parameters.

*Source: scan 4, pages 17--18.*

### 6.3 Pre-post comparison framework

Relative efficiency of the run-in test vs the standard
test: $V((\bar{C}_1 - \bar{C}_0)^*)/
V((\bar{C}_1 - \bar{C}_0)^{**}) < 1$ implies
the starred version is more efficient.

*Source: 5/13/20 and 12/29 (scan 4).*

### 6.4 ANCOVA to GLM to MMRM progression

Conceptual chain: ANCOVA -> GLM -> MMRM. The mixed model
$Y = X\beta + Zu + e$ with treatment/placebo design
matrices distinguished by the $\gamma$ column.

*Source: scan 4, page 19.*

### 6.5 General averaged change score (Wu paper)

$$d_{ik} = \frac{1}{c_i+1}\sum_{\ell=0}^{c_i}
  Y_{i,q+\ell,k}
  - \frac{1}{r_i+1}\sum_{\ell=0}^{r_i} Y_{i,\ell,k}$$

Alternative to GLS: simple $t$-test on averaged
pre-post change scores.

*Source: 7/8/21 (scan 2).*


## 7. Chronological development

| Date | Development |
|------|-------------|
| 7/27/19 | Conceptual foundations: variance components, trajectory diagrams |
| 8/25/19 | General LMM for 2-group trial, 7-step algorithm |
| 8/31/19 | Goal: reproduce Frost Eqs 11/13. Full model setup |
| 9/6/19 | Simplified model (random slope only). General variance formula |
| 9/25/19 | Change-score parameterization, normalized time |
| 10/4/19 | Woodbury for normalized time. $V^{-1} = \frac{1}{a}(I - rUU')$ |
| 10/6--7/19 | Full GLS estimator derivation (not just variance) |
| 10/13--14/19 | AR(1) sums, stratified variance, 5 correlation cases |
| 11/13/19 | Simplest run-in case: $r=1$, $k=1$ |
| 11/16/19 | General time points. Comparison of $r=1$ vs $r=0$ designs |
| 5/13/20 | Pre-post comparison framework with relative efficiency |
| 7/27/20 | Mixed model with group-specific intercepts |
| 10/15--20/20 | "Run-in notes level 3." Covariance structure development |
| 12/15--19/20 | "Long math notes." Change-score averaging, phase means |
| 3/10/21 | Test statistic, sample size formula, 5 cases formalized |
| 7/8/21 | Wu paper. General averaged change-score formula |
| 10/21 | Frost 3.1.2 by hand (consecutive changes, $r=1$, $k=1$) |
| 10/27/21 | SUR (Stanek) formulation |
| 11/11--12/21 | Frost 3.1.2: $r=2$, $k=1$ (change-from-baseline). Key result |
| 11/16/21 | General time points derivation |
| 6/23/23 | 3 time points: $(P, R, K)$ comparison framework |
| 7/20/23 | "Construct mixed model." Clean restart with Frost 3.1.1 |
| 8/27/22 | Raw-outcome Frost 3.1.1 with explicit $X$, $Z$ per group |
| 8/9/24 | Definitive Frost 3.1.1: full GLS estimator + variance |
| 8/19/24 | frost312.m finalized (Mathematica) |


## 8. What reached the manuscript

The report (`analysis/report/report.Rmd`) uses:

- The change-from-baseline parameterization (not consecutive changes)
- $R = a(I_p + \mathbf{1}_p\mathbf{1}_p')$ with Sherman-Morrison inverse
- The Woodbury identity for $\Sigma^{-1}$
- Cases 1--3: $r=0/k=2$, $r=2/k=1$, general $r/k/f$
- Numerical implementation via R code for general $(r, k, f)$
- MIRIAD parameters for AD neuroimaging examples

**Not in the manuscript** (potential extensions):

- AR(1) correlation structure and cross-covariance formulas
- Stratified variance for unequal visit patterns
- General time points (non-equally-spaced)
- The SUR formulation
- The averaged change-score ($t$-test) approach
- Random intercept contribution to variance
