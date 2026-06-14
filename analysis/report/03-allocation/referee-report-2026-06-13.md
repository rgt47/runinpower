# Referee Report: 'Optimal Allocation of Observations across Run-In, Treatment, and Common Close Phases in Longitudinal Clinical Trials'
*2026-06-13*

## 1. Summary

The manuscript (Paper 3 of a four-paper compendium) poses the design
question of how to distribute a fixed total observation budget across run-in
($J_0$), treatment ($J_1$), and common close ($J_2$) phases so as to minimise
the variance of the treatment-effect estimator $\hat\gamma$. The variance is
the $(\gamma,\gamma)$ element of $(X'\Sigma^{-1}X)^{-1}$, computed numerically
through the Woodbury machinery developed in Paper 1. The author treats the
allocation as a discrete optimisation over a finite grid, claims that the
optimum is governed primarily by the signal-to-noise ratio $\sigma_b/\sigma$,
states a $J_1$-dependent threshold above which the common close phase becomes
beneficial, and adds three extensions: non-equally-spaced timing, replicated
endpoint observations, and pattern-mixture handling of dropout following
Frost, Kenward and Fox (2008).

I read the full manuscript including every code chunk, the shared bibliography
(`analysis/references.bib`, reached through the symlink), the supporting R
source (`R/allocation.R`, `R/replicated.R`, `R/variance.R`, `R/design.R`,
`R/covariance.R`, `R/power.R`, `R/averaged.R`), the Mathematica derivations in
`analysis/scripts/*.m`, and Paper 1 (`../01-runin-power/report.Rmd`) for
cross-consistency. I re-ran the numerical claims in a clean R session that
sourced the package code directly.

## 2. Overall assessment and recommendation

**Recommendation: Reject in present form (major revision required before any
reconsideration).**

The manuscript is not yet a finished paper. It contains four unresolved `TODO`
comments at load-bearing locations (the central threshold derivation, the
general non-equal-spacing optimisation, and the entire Discussion). More
seriously, two of the three headline claims are contradicted by the paper's
own code when I re-ran it, the non-equal-spacing formula is inconsistent with
the model used in the rest of the paper, and the flagship numerical example is
degenerate (it requires one to two participants per group, so the dropout
table conveys nothing). The optimisation 'argument' is in fact an exhaustive
grid search with no analytical characterisation, which undercuts the abstract's
claim of a 'closed-form decision procedure'. Finally, the novelty claim is
overstated relative to a substantial and uncited optimal-design literature for
longitudinal mixed models.

These are correctness and substance problems, not presentation problems, so
the bar for a top-tier venue is not met at present. I detail each below with a
remedy.

## 3. Significance and novelty

The design question is real and of practical interest, particularly for the
Alzheimer's disease neuroimaging setting the compendium targets. However, the
claim in the abstract that 'the analytical optimum has not been previously
characterised' does not survive contact with the optimal-design literature.

Optimal allocation of the number and timing of measurements in linear mixed
models with a random slope is a developed field. The author cites Winkens et
al. (2005, 2007) but omits the directly relevant strand:

- Tekle, Tan and Berger and collaborators on D-optimal cohort designs for
  linear mixed-effects models, including the explicit finding that 'too many
  cohorts and repeated measurements are a waste of resources' and the result
  that the most efficient number of repeated measurements equals the number of
  cohorts plus the polynomial degree. This is the same diminishing-returns
  phenomenon the present paper rediscovers numerically.
- Fedorov and Hackl and Fedorov and Jones on optimal designs for random
  coefficient regression and multicentre trials.
- The ODmixed tool (Maastricht group) for optimal designs in heterogeneous
  longitudinal studies that explicitly incorporates dropout, which is more
  general than the pattern-mixture approximation used here.
- A 2023 paper, 'The optimal pre-post allocation for randomized clinical
  trials' (BMC Medical Research Methodology), which derives the optimal number
  of pre-treatment versus follow-up measurements and reports the general
  conclusion that repeating follow-up measurements is usually more
  advantageous than repeating pre-treatment measurements. That conclusion
  overlaps substantially with this paper's run-in-versus-replication and
  endpoint-replication sections (Sections 6.2, 6.3).

The genuinely novel element, if any, is the three-phase decomposition with a
common close phase whose treatment-effect column is constant during the close.
That is a worthwhile narrow contribution, but it is not the same as 'the
analytical optimum has not been characterised'. As written, the novelty is
overclaimed and the positioning against the optimal-design literature is
absent. Confidence: high on the existence of the prior literature (verified by
search), medium on the precise degree of overlap (I did not obtain the full
texts of every prior paper).

Searches performed (reported per instructions):
1. 'optimal design mixed effects model number timing measurements longitudinal
   D-optimal random slope Fedorov Mentre'.
2. 'Tekle Tan Berger optimal design longitudinal mixed model number repeated
   measurements allocation'.
3. 'optimal allocation pre-treatment baseline run-in measurements slope
   estimation clinical trial variance minimization'.

## 4. Major comments

### M1. The central common-close threshold claim is contradicted by the paper's own code (CORRECTNESS). VERIFIED.

The abstract and Section 5 state that the common close phase 'is beneficial
only when this ratio [$\sigma_b/\sigma$] exceeds a threshold that depends on
$J_1$; below the threshold, all observations should go to the treatment phase.'
This is the paper's flagship result.

I sourced the package and evaluated `common_close_beneficial(r=2, k, sigma2=1,
sigma_b2=ratio^2, t=1)` over a fine grid. The function returns `TRUE` (that is,
$f=1$ strictly reduces the variance over $f=0$) for **every** $k\in\{1,2,3,4\}$
and for ratios all the way down to 0.001. Direct evaluation confirms it: at
$\sigma_b/\sigma=0.1$, $k=2$, the variance falls from 0.6728 ($f=0$) to 0.4419
($f=1$). There is no positive threshold; adding a common close observation
helps at every ratio I tested.

The chunk `common_close_threshold` appears to report a threshold of 0.1 for
all $k$, but that number is an artefact: the search grid is
`seq(0.1, 5, by = 0.1)`, so 0.1 is merely the smallest value tested, and
`which(col)[1]` returns the first grid point, which is always the first. The
reported 'threshold' is the grid floor, not a real boundary. The `TODO` at the
top of Section 5.1 ('derive or compute numerically the threshold') confirms
the author never actually derived it.

Why it matters: this is the result the abstract leads with and the only
$J_1$-dependent structural claim in the paper. As stated it is false under the
paper's own model.

Remedy: either (a) prove that the common close is always weakly beneficial
under this model (which the numerics strongly suggest, since the close
observations add Fisher information about $b_i$ without cost to the treatment
contrast) and rewrite the claim accordingly, or (b) if the intended comparison
is at fixed total budget $J$ (spend a slot on close versus on treatment),
state that explicitly and re-derive. Note these are different questions:
'is $f=1$ better than $f=0$ holding $k$ fixed' (what the code computes, answer
apparently always yes) versus 'is $(J_0, J_1-1, 1)$ better than
$(J_0, J_1, 0)$' (the genuine allocation trade-off, not computed). The
manuscript conflates them.

### M2. Under fixed total budget $J$, run-in is never optimal, contradicting the framing (CORRECTNESS / SUBSTANCE). VERIFIED.

I reproduced the `optimal_grid` table and swept $\sigma_b^2$ over
$\{0.01, \ldots, 389.8\}$ and $J$ over $\{2,\ldots,10\}$. The optimal $J_0^*$
is **zero in every cell**. The optimal allocation always splits the budget
between treatment and common close (for example $J=8$, medium ratio, gives
$(0,4,4)$; the AD column gives $(0,4,4)$). Run-in observations are never chosen
when the budget is fixed.

This is a substantive and interesting finding, but it sits in direct tension
with the paper's run-in framing and with the eventual Discussion message
(currently a `TODO`) that '1-2 run-in observations ... is near-optimal'. Under
the fixed-$J$ objective the paper actually optimises, the recommendation should
be zero run-in. The reason is mechanical: with $r>0$ the model adds the third
fixed parameter $\delta$ (see `build_design`), spending a degree of freedom
without contributing to the $\gamma$ contrast, whereas run-in only helps when
it is 'free' (added on top of a fixed treatment phase, as in Paper 1's
relative-efficiency framing).

Remedy: confront this head-on. State clearly that under a fixed total-visit
budget the optimum places no observations in the run-in phase, contrast it
with the fixed-treatment-phase framing of Paper 1, and reconcile the
Discussion message. As written the paper draws the opposite practical
conclusion from what its own optimiser produces.

### M3. The 'J0 = 1 anomaly' claim is parameter-specific, not general, and the explanation is mis-stated (CORRECTNESS). VERIFIED.

Section 3.2 asserts $\Delta V_1 < \Delta V_2$ ('the first run-in observation is
less valuable than the second') 'for moderate $\sigma_b/\sigma$'. I evaluated
$\Delta V_1$ and $\Delta V_2$ across ratios:

| $\sigma_b/\sigma$ | $\Delta V_1$ | $\Delta V_2$ | $\Delta V_1 < \Delta V_2$? |
|---|---|---|---|
| 0.10 | 0.2605 | 0.0867 | no |
| 0.32 | 0.1633 | 0.0080 | no |
| 0.50 | 0.0600 | 0.0259 | no |
| 1.00 | 0.1765 | 0.6667 | yes |
| 2.00 | 4.2000 | 2.1718 | no |
| 5.00 | 44.74 | 3.444 | no |

The inequality holds at exactly one of the tested ratios ($\sigma_b/\sigma=1$),
which is the value used in the displayed table ($a=b=t=1$). It fails at the
paper's own 'low' ratio (0.32) and at high ratios. So the 'anomaly' is an
artefact of the single parameter setting in the table, not a 'moderate-ratio'
phenomenon. The phrase 'for moderate $\sigma_b/\sigma$' is unsupported.

The degrees-of-freedom explanation ('the first run-in observation forces the
model from 2 to 3 parameters ... the second does not increase the parameter
count and provides pure information gain') is a reasonable intuition for why
$\Delta V_1$ can be small, but it does not establish $\Delta V_1 < \Delta V_2$
in general, and the numerics show it does not hold in general. The analogy to
Frison and Pocock (1992) is loose.

Remedy: restrict the claim to the specific parameter regime where it holds,
state that regime quantitatively, and either prove the boundary or present the
table above. Do not generalise from $n=1$ parameter setting.

### M4. The non-equal-spacing variance formula is inconsistent with the paper's model and is not verified (CORRECTNESS). VERIFIED (mismatch); SUSPECTED (error).

Section 6.1 displays
$\Var(\hat\gamma) = 3\sigma^2/(1 - t_1 t_2) + 2\sigma_b^2$
'from the notebook derivations' under $\sum t_i^2 = 1$, $\sum t_i = 0$. Three
problems:

1. Model mismatch. The constraints $\sum t_i = 0$ (centred times) and the
   diagonal residual $R = \sigma^2 I$ come from
   `analysis/scripts/raw_outcome_model.m`, which uses raw outcomes with a
   random slope and no random intercept. The rest of the paper (and `R/`)
   uses change scores from a fixed baseline with $R = \sigma^2(I + \mathbf 1
   \mathbf 1')$ and design times $t, 2t, \ldots$ that are strictly positive
   (so $\sum t_i = 0$ is impossible for change scores). The spacing section is
   therefore in a different model from everything else and is never reconciled.

2. The formula does not match my from-scratch GLS computation in the
   change-score model. For two change scores at times $(u_1, u_2)$,
   `var_gamma_matrix(0,2,0,1,1,1)` gives 3.0, while a direct GLS evaluation at
   $(1,2)$ gives 1.5; the displayed formula gives neither.

3. In the raw-outcome model it claims to come from, the formula is only
   approximately correct. At the symmetric optimum $t_1 = -t_2 = 1/\sqrt 2$ I
   get a two-group slope variance of exactly 4.0, matching the formula. But at
   $(0.6, -0.8)$ I get 4.0408 while the formula gives 4.0270. The discrepancy
   suggests the closed form is either an approximation or carries a small
   algebra error.

The accompanying statement 'minimized when $t_1 t_2$ is minimized (most
negative)' has the right direction (minimising variance means maximising
$1 - t_1 t_2$) but is argued from a formula I could not verify.

Remedy: state the model for this section explicitly, re-derive the closed form
in that model, verify it numerically against a direct GLS computation at
several non-symmetric configurations, and reconcile (or clearly separate) it
from the change-score formulation used elsewhere. As it stands the only fully
general claim in the non-equal-spacing section is followed immediately by a
`TODO` admitting the general case was not done.

### M5. The flagship AD numerical example is degenerate (SUBSTANCE). VERIFIED.

With the MIRIAD parameters ($\sigma^2 = 0.0025$, $\sigma_b^2 = 0.9745$,
$t = 1$, $\Delta = 0.25$), `sample_size(2,2,1,...)` returns $N = 1$ per group,
and the partial-pattern stratum returns $N = 2$. Consequently the entire
dropout table (Section 8) reads $N = 1$ at 0 percent dropout and $N = 2$ at
10, 20 and 30 percent dropout. A required sample size of one or two
participants per arm is not a usable illustration; it tells the reader nothing
about how dropout interacts with allocation, which is the table's stated
purpose. The cause is the very small $\sigma^2$ relative to $\Delta^2$
($(z_\alpha + z_\beta)^2 v_1 / \Delta^2 \approx 10.5 \times 0.0049 / 0.0625$).

Remedy: choose $\Delta$, $\sigma^2$ and $\sigma_b^2$ that yield realistic
sample sizes (tens to hundreds per arm), or rescale the effect size to the
units actually used in MIRIAD-type trials. Re-run the dropout table so that it
exhibits non-trivial variation. Also reconcile with Paper 1, which presumably
uses the same parameters; if Paper 1's $N$ is also degenerate, that is a
compendium-wide problem.

### M6. The pattern-mixture dropout treatment is mis-specified at the stratum level (CORRECTNESS). INSPECTED.

In the `dropout` chunk, the 'partial' stratum is computed as
`sample_size(r_opt, k_opt, 0, ...)`, that is, the dropout stratum keeps the
full $k=2$ treatment observations and merely loses the common close. That is
not what dropout does: a participant who drops out loses the *later* visits,
so the realistic partial pattern has fewer treatment observations (smaller
$k$, or a truncated visit set), not the same $k$ with $f=0$. As coded, the
'dropout' stratum is simply the $(2,2,0)$ design, which is a design choice, not
a dropout pattern. The pattern-mixture combination formula
$N = 1 / \sum_s p_s / N_s$ is then applied to strata that do not correspond to
the dropout process described in the text.

Why it matters: the section claims to operationalise the Frost-Kenward-Fox
pattern-mixture approach, but the strata are not dropout patterns, so the
numbers do not represent dropout-adjusted sample sizes.

Remedy: define each completion stratum by the visits actually observed under
monotone dropout (for example, dropout after visit $j$ yields the design with
the first $j$ post-baseline observations), compute $N_s$ for each such
truncated design, and weight by the probabilities of those patterns. The
manuscript already cites Hu, Mackey and Thomas (2021) as the more rigorous
alternative; if the pattern-mixture approximation is retained it must at least
use genuine dropout strata.

### M7. The 'optimisation' is an exhaustive grid search presented as a closed-form procedure (EXPOSITION / SUBSTANCE). VERIFIED.

`optimal_allocation` enumerates all $(r,k,f)$ with $r+k+f=p$ and returns the
argmin. There are no first-order conditions, no KKT analysis, no monotonicity
lemmas, and no characterisation of the optimum as a function of
$\sigma_b/\sigma$ beyond reading values off a table. The abstract nonetheless
advertises 'a closed-form decision procedure' and 'closed-form tractability'.
A brute-force search over a small grid is neither closed form nor a
characterisation; it is a computation. For a methods journal at this tier, the
reader expects either an analytical optimum (sign of a derivative, a proven
threshold, a monotone comparative static) or, at minimum, an explicit argument
that the discrete objective is unimodal so the grid search is guaranteed
correct. Neither is present.

Remedy: derive the comparative statics. The variance is a ratio of low-degree
polynomials in the design counts; with the Woodbury structure it should be
possible to sign $\partial V / \partial f$ and $\partial V / \partial r$ and to
prove the claimed dependence on $\sigma_b/\sigma$ analytically. Failing that,
reframe the contribution honestly as a numerical exploration and drop
'closed-form' from the abstract and conclusions.

## 5. Minor comments

- The manuscript contains four `TODO` comments at substantive locations:
  Section 5.1 (the threshold derivation, see M1), Section 6.2 (the general
  non-equal-spacing optimisation), and the entire Discussion (Section 9). A
  paper with an empty Discussion and an undone central derivation is not
  submittable. INSPECTED.

- The abstract states non-equal spacing improves efficiency '5-15 percent ...
  with the largest gains at low signal-to-noise ratios', but no chunk in the
  manuscript computes this range; the only spacing content is the single
  two-point formula in 6.1 and a `TODO`. The 5-15 percent figure is
  unsupported by anything in the document. SUSPECTED (claim without backing
  computation).

- The abstract claims replication has 'diminishing returns above 3-4
  replicates'. The replication table (`replication_table`) reduces variance by
  the algebraic factor governing $\sigma^2/n_{\text{rep}}$; the largest
  fractional reduction I reproduced is 0.31 at $n_{\text{rep}}=5$,
  $\sigma_b/\sigma=0.5$. 'Diminishing returns' is correct qualitatively, but
  the specific '3-4' cutoff is not derived. Tie the claim to a computed
  quantity. VERIFIED (table reproduces); SUSPECTED (the 3-4 cutoff).

- `var_gamma_replicated` contains dead code: lines that build `R` via
  `diag(sigma2_vec) + sigma2/n_rep_pre * matrix(1,p,p)` and then immediately
  overwrite the diagonal, followed by a separate `R_correct` loop that is the
  matrix actually used. The first `R` is never used. I verified the function
  reduces correctly to `var_gamma_matrix` at $n_{\text{rep}}=1$ for several
  configurations (max difference $0$), so the result is right, but the dead
  code should be removed for the reproducibility artefact to be clean.
  VERIFIED.

- In `var_gamma_replicated`, the `if (f > 0) ... else ...` branch (lines 46-50)
  assigns the same value `sigma2 / n_rep_post` to `sigma2_vec[p]` in both
  branches, so the conditional is vacuous. More importantly, when $f>0$ the
  last visit is a *common close* visit, not the treatment endpoint, so
  replicating index $p$ does not replicate 'the final treatment visit' the
  text describes. The replication target is mislabelled whenever $f>0$.
  INSPECTED.

- Reproducibility section inconsistencies: the prose refers to
  'analysis/paper{1,2,4}-*/' and 'analysis/paper3-allocation/', but the actual
  directories are `analysis/report/0{1,2,3,4}-*`. The Data Availability
  statement lists `design.R` among the optimal-allocation routines, but
  `design.R` builds design matrices; the allocation routines are in
  `allocation.R` and `replicated.R`. The 'set.seed(20260418)' note refers to
  Monte Carlo benchmarking that does not appear in this manuscript. INSPECTED.

- Geometry `right=5cm` produces an unusually wide right margin; likely a
  copy-paste artefact from a template. Cosmetic.

- The compound-symmetry correlation in Eq. (between lines 206-209),
  $\Corr(C_{-J_0}, C_l) = \sigma_b^2 t^2 / (2\sigma^2 + \sigma_b^2 t^2)$, is
  stated as 'constant in $J_0$'. It is also constant in $l$, which makes the
  subsequent sentence ('determines the information gain') hard to interpret,
  since a constant correlation cannot by itself explain a $J_0$-varying
  marginal value. The author in fact concedes the diminishing returns come
  from the information-matrix structure, not the correlation, which raises the
  question of why the correlation display is included at all. Clarify.
  INSPECTED.

- Notation versus Paper 1: Paper 1 uses $J_0, J_1, J_2$ in prose and $r, k, f$
  in code; Paper 3 does the same. This is internally consistent. The core
  variance for the conventional design, $(\sigma^2 + 2\sigma_b^2 t^2)/t^2$
  (Paper 1 Eq. var_frost, and `var_gamma_frost`), is consistent with the
  Mathematica Part 1 result $(a + 2 b t^2)/t^2$. The shared variance engine
  `var_gamma_matrix` is identical across papers. Cross-paper consistency on the
  variance formula and notation is satisfactory. VERIFIED.

## 6. Missing or questionable references

| Citation | Issue | Action |
|---|---|---|
| Tekle, Tan, Berger (D-optimal cohort designs for LMM; 'too many cohorts and repeated measurements are a waste of resources') | Directly relevant prior work on optimal number/timing of measurements in random-slope LMMs; not cited | Add and position against |
| Fedorov and Hackl (1997); Fedorov and Jones (2005) | Foundational optimal-design results for random coefficient regression and trials; not cited | Add |
| 'The optimal pre-post allocation for randomized clinical trials' (BMC Med Res Methodol, 2023) | Directly overlaps Sections 6.2-6.3 (pre vs post measurement allocation, replication) | Add and reconcile claims |
| Mentre, Mallet, Baccar (population optimal design) | Standard reference for optimal design in mixed models | Add |
| Atkinson and Donev (Optimum Experimental Designs) | Standard textbook for the optimisation framing the paper claims to introduce | Add |
| ODmixed (optimal designs for longitudinal studies with dropout) | More general than the pattern-mixture approximation in Section 8 | Add and compare |
| `hu2021` (bib) | Defined three times (lines 255, 505, 701) with two different papers sharing one key; BibTeX silently keeps one, so the Section 8 citation is ambiguous | Split into distinct keys (e.g. hu2021bimj, hu2021ijb) |
| `nash2021`, `iddi2022` | Cited correctly as implementations that do not optimise allocation; fine, but the Stata `slopepower` and `longpower` framing should note they do address number of measurements indirectly | Minor clarification |
| `goulet2019` | Cited for '20-50 replications' diminishing returns; this is a psychology-methods paper on replicate measurements, a reasonable but not domain-matched citation; verify the quoted range | Verify |

Bib hygiene: the triple `hu2021` key is the only hard error and must be fixed.
The file is otherwise a shared 60-entry collection; many entries are unused in
this manuscript (acceptable for a shared bib, but consider a per-paper subset).

## 7. Suggestions

1. Resolve the question in M1 analytically: prove whether the common close is
   always weakly beneficial at fixed $k$, then separately characterise the
   genuine fixed-$J$ trade-off. This single result, done rigorously, could
   anchor the paper.
2. Lead with the M2 finding (zero run-in under fixed budget) as the honest
   headline; it is counterintuitive and useful, and it differentiates the
   fixed-budget question from Paper 1's fixed-treatment-phase framing.
3. Replace the grid search with comparative statics (sign of
   $\partial V/\partial f$, $\partial V/\partial r$) so the paper delivers the
   characterisation its abstract promises.
4. Fix the AD example parameters so sample sizes are realistic, and rebuild the
   dropout table on genuine monotone-dropout strata (M5, M6).
5. Position explicitly against the optimal-design literature in Section 1;
   state precisely what the three-phase / common-close structure adds beyond
   Winkens, Tekle-Tan-Berger, and the 2023 pre-post allocation paper.
6. Remove all `TODO`s; write the Discussion; either complete or remove the
   general non-equal-spacing section rather than leaving a single special case
   plus a placeholder.

## 8. Scope of review

I verified the numerical claims by sourcing the package R code
(`R/*.R`) in a clean R 4.5.3 session and re-running the manuscript's
computations plus targeted sensitivity sweeps; results labelled VERIFIED were
reproduced directly. I inspected the Mathematica derivations
(`analysis/scripts/woodbury_variance.m`, `rcrm_raw_outcome.m`,
`raw_outcome_model.m`) and read the source code paths flagged above; claims
labelled INSPECTED rest on reading rather than execution. I did not run the
full `rmarkdown` render (the build is container-based and the `tools/`
stage-render harness was not exercised), so I cannot certify that the typeset
PDF matches the chunk outputs I reproduced, though the chunk code is what I
ran. I confirmed the existence and relevance of the prior optimal-design
literature by web search (three queries, reported in Section 3) but did not
obtain every full text, so the precise degree of overlap is assessed at medium
confidence. I did not independently re-derive the Woodbury variance formula
from Paper 1 (that is the subject of the Paper 1 review); I took it as given
and checked only its consistent use here.

---
*Rendered on 2026-06-13 at 13:00 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/analysis/report/03-allocation/referee-report-2026-06-13.md*
