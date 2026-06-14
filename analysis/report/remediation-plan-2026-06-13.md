# Remediation Plan: Run-In / Common-Close Compendium
*2026-06-13 18:22 PDT*

Concrete, anchored fixes for all four manuscripts, derived from the
four individual referee reports and the pooled report in this
directory. Each item gives the location, the change, and a
verification step. Work top to bottom: Phase 0 changes propagate
into every paper and must be done first, because they alter the
numbers that Phases 1 to 4 then report.

Legend for confidence on the underlying finding: [V] verified by
re-running the paper's own code, [I] inspected, [S] suspected.

---

## Phase 0. Shared infrastructure (do once, affects all four)

### 0.1 Fix the MIRIAD / Alzheimer's parameter set [V, highest value]

Problem. All four papers set residual `sigma^2 = 0.0025` (Frost
Table II 'additional measurement error') with `sigma_b^2 = 0.9745`
and silently drop `sigma_u^2 = 0.3338` (within-subject between-visit
variance), which is the component that plays the residual role in a
change-score analysis. With `Delta = 0.25` on this scale the
examples collapse to N = 1 to 3 per arm everywhere.

Fix.

1. Create one canonical parameter block and source it from all four
   papers (for example `analysis/scripts/miriad_params.R` or a
   shared chunk):
   - `sigma_b2 <- 0.9745`  (between-subject slope variance)
   - `sigma_u2 <- 0.3338`  (within-subject between-visit)
   - `sigma_meas2 <- 0.0025`  (additional measurement error)
   - change-score residual `sigma2 <- sigma_u2 + sigma_meas2`
     (= 0.3363), with an explicit derivation of how `sigma_u^2`
     enters `R` for the change-score model (Paper 1 M2).
2. Re-derive `R` if the model is meant to omit the visit-level
   random term; if so, justify the omission and show its effect on
   `R = sigma^2(I + 11')` (Paper 1 M2 remedy).
3. Re-run every numerical example, table, and figure in all four
   papers against the corrected block.

Verification. Paper 1 should reproduce Frost (2008) Table III as an
external check: roughly 175 per arm (12-month, no run-in) and 122
(24-month). With residual 0.3363 the Paper 1 reviewer obtained
n = 385 -> 145 (62% reduction), a sane range. Confirm Papers 2 to 4
now yield tens-to-hundreds per arm, not 1 to 3.

Downstream consequence. Every reported percentage changes. Do not
edit the abstract numbers until this is done (see 1.3, 2.4, 3.5,
4.5).

### 0.2 Repair the shared `references.bib` [V]

Problem. One `analysis/references.bib` is symlinked into all four
papers and is a Zotero dump with hard errors.

Fix.

1. `hu2021` is defined two-to-three times (lines 255, 505, 701)
   for two different real papers under one key. BibTeX silently
   keeps one, so at least one paper mis-cites. Split into distinct
   keys and cite each for its result:
   - `hu2021bimj`: Hu N, Mackey HM, Thomas R. Biom J.
     2021;63(4):806-824 (monotone-missing).
   - `hu2021ijb`: Hu N, Mackey H, Thomas R. Int J Biostat.
     2021;18(1):173-182 (power formula).
2. De-duplicate entry pairs for the same work: keep one of
   `frost2008` / `frostOptimizingDesignClinical2008`,
   `wang2019` / `wangTwoPeriodLinear2019`,
   `laird1982` / `lairdlRandomEffectsModelsLongitudinal2023`.
   The `laird...2023` entry has corrupted authors ('Lairdl, Nan M',
   'Warel, James H') and a wrong year; delete it and cite the
   correct `laird1982` (Biometrics 1982).
3. Standardise the citation-key scheme across the compendium (Paper 1
   uses `@frostOptimizingDesignClinical2008`; Paper 2 uses
   `@frost2008`). Pick one key per work and apply everywhere.
4. Verify DOIs: `frost2008` second block (line 483) has wrong issue
   number = {18}; the correct entry has number = {19}.
   `winkens2005` DOI `10.1002/sim.2370` may be wrong (web record
   shows `10.1002/sim.2385`); verify.
5. Prune the ~50 uncited entries, or maintain a per-paper subset at
   submission.

Verification. `grep '@' references.bib | sort` shows no duplicate
keys; a test render of each paper produces no BibTeX 'repeated
entry' warning and all `\cite` keys resolve.

### 0.3 Standardise notation across the compendium [I]

Problem. Prose uses `J0/J1/J2`; all code and the R package use
`r/k/f`. Symbol collisions compound it: `a,b` denote both the random
effects `a_i,b_i` and the variances `sigma^2,sigma_b^2` (Paper 1
line 360); Paper 2 redefines `c = rho` (line 271) whereas Paper 1
and the Mathematica sources use `c = sigma^2 - sigma_b^2 t^2`.

Fix.

1. Choose one design-count convention and use it in prose, equations,
   and code. If `J0/J1/J2` stays in prose, add a single explicit
   mapping line (`J0 = r, J1 = k, J2 = f`) and rename code, or
   vice versa (Paper 1 M4, Paper 4 m1).
2. Rename the variance shorthands so they do not collide with the
   random effects `a_i,b_i`; the Mathematica `c,d` reparameterisation
   may live in an appendix (Paper 1 M4).
3. Give the AR(1) autocorrelation in Paper 2 a distinct symbol (not
   `c`), reserving `c` for the Paper 1 reparameterisation (Paper 2
   minor / pooled 2.3).
4. Add a single symbol table, shared across the four papers.

### 0.4 Fix reproducibility metadata across all papers [V/I]

1. Repository paths are wrong project-wide: papers cite
   `analysis/paperN-*/`; the actual layout is `analysis/report/0N-*/`.
   Correct the Data Availability and Reproducibility sections in all
   four (Paper 2 minor, Paper 3 minor, Paper 4 m6).
2. Derivation-source provenance: CLAUDE.md and the manuscripts name
   `frost311.m` / `frost312.m`, which do not exist. Actual files:
   `woodbury_variance.m`, `blup_derivation.m`, `raw_outcome_model.m`,
   `rcrm_raw_outcome.m`. The scripts are Mathematica/Wolfram but are
   called 'MATLAB' (Paper 1 lines 1207, 1233, 1237). Correct names
   and language; update CLAUDE.md (Paper 1 M8).
3. Remove fabricated or absent Monte Carlo claims wherever no
   simulation exists (Paper 1 lines 1239-1243 `set.seed(20260418)`;
   Paper 2 lines 268-270; Paper 4 M3, see 4.3). Either implement the
   simulation or delete the claim.
4. Reconcile testing-harness references: Paper 4 cites `tinytest`;
   the project uses `testthat` (Paper 4 m5).

---

## Phase 1. Paper 1 (01-runin-power) -- the anchor

Recommendation: major revision. Verified algebra is sound; fix the
two false/overstated claims, complete and validate the common-close
contribution, then bring this paper to submittable as the anchor.

### 1.1 Correct the false J0 = 1 theorem [V] (Section 6, lines 934-944)

The claim that `Var|J0=1` can exceed `Var|J0=0` is false for this
estimator. Exhaustive grid (k in 1..6, sigma_b^2 in [1e-4,1e3], t in
{0.25,0.5,1,2,4}): J0=1 is never worse, and exactly equal at
sigma_b^2=0.5, t=1, k=2. Replace with the correct 'never worse,
sometimes no gain' statement plus the boundary condition for zero
gain. Drop or re-scope the Frison-Pocock (1997) citation (different
estimator and parameterisation). Remove the 'plan zero or two,
never one' recommendation.

### 1.2 Rebuild Section 5.5 AD application [V] (lines 1093-1156)

Driven by 0.1. After fixing the residual, re-run all of 5.5 and the
relative-efficiency figure. Reconcile against Frost Table III. State
explicitly how `sigma_u^2` maps into the change-score `R`.

### 1.3 Correct the abstract numbers [V] (lines 29-44)

The '30-45%' figure is Frost's MIRIAD headline, not this paper's
result (the paper's own run-in table gives 14/28/76%). Replace with
figures from the corrected analysis (1.2) and state the parameter
regime. Delete the unverifiable 'simplest available closed-form'
superlative, which also conflicts with the paper's own line 808
('does not admit a simple closed-form expression').

### 1.4 Deliver and validate the common-close contribution [I]

This is the paper's one clearly novel element and is currently a
numerical recipe only.

1. Derive a genuine closed form for at least J2 = 1 and J2 = 2
   (Section 7 suggestion 1).
2. Add a direct (non-Woodbury) GLS evaluation
   `Var = [(X' Sigma^{-1} X)^{-1}]_{gamma,gamma}` with Sigma
   assembled directly, plus a small Monte Carlo, for several
   (J0,J1,J2) with J2 > 0. No current check exercises J2 > 0 (M6);
   the two existing checkpoints both have J2 = 0.
3. Derive the common-close change-score mean from the generative
   model rather than asserting `g_i = 0`; state the carryover /
   washout / level-maintenance assumptions and the estimand; contrast
   with delayed-start designs (Liu-Seifert 2015; Wang 2019) (M5).

### 1.5 Supply the missing r = 2 derivation trail [V/S] (M7, M8)

`woodbury_variance.m` Part 2 derives the r = 1 case, not the
headline r = 2 case; `var_r2k1` (Eq., lines 666-671) is asserted with
a cofactor sketch the reviewer could not reproduce symbolically
(though it matches the matrix routine to 3.6e-15 across 12 sets).
Supply the symbolic r = 2 derivation or cite the exact script.
Either prove every entry of Eq. XRX_general (lines 760-773) or
downgrade it to 'the implementation computes' (M7).

### 1.6 Paper 1 minor fixes

- Strip hard-coded absolute paths and the user-specific render hook
  (line 56, lines 1245-1250).
- Fix the change-score time list for J=3 (line 474 lists four times
  for three change scores; the j=0 baseline is differenced away).
- Delete the malformed `Z = ...` expression (lines 487-492).
- Rewrite Eq. H (lines 747-750) using `eta_J` to remove circularity.
- Figure axis labels use literal `J_0`, `J_2` underscores in base R;
  use `expression()` or ggplot (line 1065).
- Add missing citations: Nash et al. 2021 slopepower (+ 2024
  erratum); elevate Iddi-Donohue 2022 longpower to a numerical
  benchmark; cite Galbraith-Marschner 2002 if its guidance is used.

---

## Phase 2. Paper 2 (02-correlation)

Recommendation: reject in present form / major revision. Complete the
draft, then fix the formula error and the two false headline claims.

### 2.1 Fix the moment-based variance-of-mean formula [V] (M2, lines 291-308)

The lag sum uses `c^{i-1}` where `c^{i}` is required; off by a factor
`c` per term (n=2, rho=0.3: returns 1.000 vs correct 0.650). Replace
with `Var(Xbar_n) = (sigma^2/n^2)(n + 2 T_n)`,
`T_n = sum_{l=1}^{n-1}(n-l)c^l = c[n(1-c)-(1-c^n)]/(1-c)^2`. The
cross-covariance term `ar1_cov_phase_means` is correct; only the
variance term and the inherited `ar1_var_change` are wrong. Re-run
Section 4 and all dependent output. Verify against
`1' Sigma 1 / n^2`.

### 2.2 Resolve the two-route equivalence claim [V] (M1, abstract lines 22-25)

The abstract claims the Woodbury and moment routes 'yield identical
results'; the only numerical comparison (chunk `moment_based`, lines
295-326) prints 1.4375 vs 1.6662, and Section 5 (lines 329-353) says
they 'compute different quantities'. They target different estimands
(slope-difference gamma vs a phase-mean difference). Either (a)
define a common estimand (set sigma_b^2 = 0 and prove algebraic
equality), or (b) withdraw the equivalence claim and present the
moment route as a distinct estimator. Note (a) needs 2.3 first.

### 2.3 Fix `var_gamma_ar1` at sigma_b^2 = 0 [V] (M3, lines forming H)

`solve(solve(G))` is singular at sigma_b^2 = 0, blocking the only
valid equivalence check. Use the equivalent
`H = G - G Z'(R + ZGZ')^{-1} Z G`, or guard the limit analytically.

### 2.4 Correct the 5-25% sample-size claim [V] (M4, abstract lines 32-36)

Running the paper's own `comparison` chunk (MIRIAD, rho=0.5): AR(1)
sizes are lower than or equal to CS, not 5-25% higher, and the
deviation grows with rho, so 'largest at small rho' is backwards.
Recompute and re-state. Decide the cross-structure matching
convention (fixed nominal rho vs common lag-1 correlation vs common
total information) and characterise the sign of the effect correctly.
The implausible N = 1 to 3 is the 0.1 units problem.

### 2.5 Complete the unfinished sections [I] (M7, M8)

Three sections carry `TODO`s, the AR(1) change-score correlation
derivation reads '(to be derived)' (lines 172-175), and the
Discussion is a single `TODO` (lines 475-476). Either derive the
five-case taxonomy boundaries (the (rho, sigma_b/sigma, J0) region
where CS and AR(1) agree to within a tolerance) or drop the taxonomy
framing (M8). Write the Discussion.

### 2.6 Paper 2 minor / positioning

- Fix the change-score definition inconsistency (Section 2.2 declares
  `C_j = Y_j - Y_0` but the formulas are correct only for
  `C_j = Y_j - Y_{j-1}`; AR(1) code uses change-from-baseline)
  (M6, lines 159-170).
- Correct the false novelty premise: `longpower::power.mmrm.ar1`
  (Lu-Luo-Chen 2008) is exactly the analytic AR(1) formula the intro
  (lines 114-127) says does not exist. Reframe around the narrower
  defensible claim (AR(1) + run-in/common-close + random slope).
- Cite and position against Frison-Pocock 1992 (in bib, uncited),
  Munoz et al. 1992 (missing), Winkens et al. 2007 (engage, not in
  passing).
- Fix the 'CS is a limit of AR(1)' conflation (M5, lines 124-127):
  the model collapses to CS at rho -> 0 (no serial correlation), not
  a compound-symmetric limit.
- State `var_gamma_ar1` correctness (matches block-diagonal GLS to
  machine precision) as the paper's one solid result.
- Fix the `geometry: right=5cm` template artefact.

---

## Phase 3. Paper 3 (03-allocation)

Recommendation: reject in present form / major revision. Resolve the
threshold question analytically, confront the fixed-budget result,
and reposition against the optimal-design literature.

### 3.1 Resolve the common-close threshold [V] (M1, Section 5, abstract)

`common_close_beneficial()` returns TRUE for all k=1..4 down to
ratio 0.001; the reported 'threshold 0.1' is the grid floor
(`seq(0.1,5,by=0.1)` with `which(col)[1]`). No positive threshold
exists. Either (a) prove the common close is always weakly beneficial
at fixed k and rewrite the claim, or (b) if the intended comparison
is at fixed total budget (a close slot vs a treatment slot), state
and re-derive that. The paper currently conflates
'f=1 vs f=0 at fixed k' (always yes) with
'(J0,J1-1,1) vs (J0,J1,0)' (the genuine trade-off, not computed).

### 3.2 Confront the fixed-budget result [V] (M2, chunk `optimal_grid`)

Under fixed total J, optimal J0* = 0 in every cell; the budget always
splits between treatment and common close. This contradicts the
run-in framing and the intended Discussion message ('1-2 run-in
near-optimal'). State plainly that under a fixed-visit budget the
optimum places no run-in observations, and contrast with Paper 1's
fixed-treatment-phase framing. Consider leading with this as the
honest headline.

### 3.3 Restrict the J0 = 1 anomaly claim [V] (M3, Section 3.2)

`DeltaV_1 < DeltaV_2` holds only at sigma_b/sigma = 1 (the table's
single setting); it fails at the paper's own 'low' ratio 0.32 and at
high ratios. Restrict the claim to the regime where it holds, state
that regime quantitatively, and either prove the boundary or present
the full table. Drop 'for moderate sigma_b/sigma'.

### 3.4 Fix the non-equal-spacing formula [V mismatch, S error] (M4, Section 6.1)

`3 sig^2/(1 - t1 t2) + 2 sigb^2` comes from `raw_outcome_model.m`
(raw outcomes, random-slope-only, centred times, R = sig^2 I) while
the rest of the paper uses change scores with R = sig^2(I + 11') and
strictly positive times (so sum t_i = 0 is impossible). It does not
match direct GLS and is only approximate even in its own model
(4.0408 vs 4.0270 at (0.6,-0.8)). State the model for this section,
re-derive in it, verify numerically against direct GLS at
non-symmetric configurations, and reconcile or clearly separate it
from the change-score formulation. Then either complete the general
optimisation or remove the section (do not leave the `TODO`).

### 3.5 Fix the degenerate example and dropout strata [V/I] (M5, M6)

Driven by 0.1. After the parameter fix, rebuild the dropout table so
it varies. Redefine each completion stratum by the visits actually
observed under monotone dropout (dropout after visit j -> design with
the first j post-baseline observations), compute N_s per truncated
design, and weight by pattern probabilities. The current 'partial'
stratum keeps full k=2 and only drops f, which is a design choice,
not a dropout pattern (M6).

### 3.6 Reframe the 'optimisation' honestly [V] (M7)

`optimal_allocation` is an exhaustive grid search; the abstract calls
it a 'closed-form decision procedure'. Either derive comparative
statics (sign of dV/df, dV/dr; a proven threshold; a monotonicity /
unimodality argument) or drop 'closed-form' and present it as a
numerical exploration.

### 3.7 Paper 3 minor / positioning

- Remove the unsupported abstract figures with no backing
  computation: non-equal spacing '5-15%' (no chunk computes it) and
  replication 'diminishing returns above 3-4 replicates' (the 3-4
  cutoff is not derived). Tie each claim to a computed quantity.
- Clean `var_gamma_replicated`: remove the dead `R` construction and
  the vacuous `if (f>0) ... else ...` branch (lines 46-50); fix the
  mislabelled replication target when f > 0 (replicating index p is a
  common-close visit, not the treatment endpoint).
- Position against the optimal-design literature: Tekle-Tan-Berger
  (D-optimal cohort designs), Fedorov-Hackl / Fedorov-Jones, ODmixed
  (optimal design with dropout), the 2023 BMC MRM optimal pre-post
  allocation paper, Mentre et al., Atkinson-Donev. Drop 'the
  analytical optimum has not been previously characterised'.
- Clarify the constant-correlation display (lines 206-209) that
  cannot by itself explain a J0-varying marginal value.

---

## Phase 4. Paper 4 (04-gls-vs-ttest)

Recommendation: reject in present form. The ordering result is the
only surviving headline; the other two are false and the simulation
does not exist.

### 4.1 Correct the central ANCOVA '<5% loss' claim [V] (M1, abstract lines 35-39, chunk `re_table`)

Re-running the paper's own grid: RE(GLS/ANCOVA) is mostly 0.1-0.5
(ANCOVA variance 2-10x GLS); the >=0.95 condition holds in only 5.2%
of cells, all with J1=1. Root cause: `var_gamma_ancova` collapses the
k post-visits into one mean, discarding the slope information GLS
uses. The formula is correctly implemented (MC ratio 1.007); the
verbal claim is wrong. Either delete the '<5%' claim (and the
'ANCOVA is the recommended default' conclusion), or replace
ANCOVA-on-means with an ANCOVA that retains all post visits
(per-visit response with run-in covariate, or cLDA) for which a
near-GLS claim might hold.

### 4.2 Correct or retract the robustness claim [V] (M2, chunks `misspec_G`, `misspec_R`)

The figures show the opposite of the claim: misspecified-GLS variance
peaks at 5.3 (misspec_G) and 0.74-2.16 (misspec_R) while averaged is
~20 throughout; GLS keeps a ~10x margin even when misspecified.
Either retract, or construct a misspecification regime where a simpler
estimator actually overtakes GLS (compare misspecified-GLS variance
against a correctly analysed simpler estimator including t-test
small-sample behaviour) and report the crossover.

### 4.3 Implement or remove the Monte Carlo study [V] (M3, M4, Section 'Simulation study')

The text describes a 5000-replicate factorial MC with empirical
variance, power, and coverage; the `simulation` chunk runs ONE trial
under `set.seed(42)`. The Reproducibility section cites
`set.seed(20260418)` and per-replicate `+ rep_idx` seeding that
appear nowhere in the source, and names `nlme` as the GLS fitter
though `nlme` is never called. Either implement the full factorial
simulation with Monte Carlo standard errors and a
formula-vs-simulation validation table, or delete the simulation
section and every methods/abstract sentence that references it. This
is a fabricated-methods correction and is mandatory before any
submission.

### 4.4 State and prove the ordering result [V] (M6)

GLS <= ANCOVA <= averaged holds with zero violations across 413 grid
points, but the proof is only gestured at and the section that should
contain it is an empty `TODO` (line 405). Supply a proposition: both
simpler estimators are linear unbiased estimators of gamma under the
random-slopes model, so each has variance at least the GLS/BLUE; and
ANCOVA <= averaged because ANCOVA optimises over the run-in
coefficient phi while averaged fixes it. Characterise the boundary
analytically.

### 4.5 Fix the degenerate MIRIAD example [V] (M5, chunk `miriad_comparison`)

Driven by 0.1. At sigma_b/sigma ~ 19.7 the table gives GLS N = 1
alongside averaged N up to 8193, uncommented, and contradicts 4.1.
After the parameter fix, justify components and Delta on the original
scale and re-run.

### 4.6 Complete or remove unfinished sections, add cLDA [I] (M7, M8)

The SUR subsection is a stub + `TODO`, 'When to use which estimator'
is a `TODO`, and the Discussion is a comment block (~1/3 of intended
content absent). Complete or remove each. Add the cLDA estimator
(Liang-Zeger 2000; Lu 2010; Lu-Mehrotra-Liu 2009 for sample size):
it sits between averaged and GLS, exploits the covariance without
specifying G, and is the natural robust comparator that would deliver
the near-GLS-with-robustness story the abstract wants. It is arguably
the point of the comparison.

### 4.7 Paper 4 minor / positioning

- Cite Frison-Pocock 1992/1997 (in bib, NOT cited in this paper
  despite being the closest prior work) and Liang-Zeger 1986 at the
  sandwich construction. Verify the crager1987 attribution.
- Tighten 'MVUE' to BLUE (abstract, m7).
- Reword the misspec figure captions: the averaged estimator's
  variance does change with the true R; it is the estimator
  definition that does not depend on the assumed model (m8).
- Fix the subject-varying limits in Eq. (avg) vs the fixed-design
  development (m2); state the population-ANCOVA (phi known) caveat
  and the feasible-estimator inflation (m3).
- Remove the leftover single-trial `set.seed(42)` debugging output
  (m9).

---

## Suggested execution order

1. Phase 0 in full (parameters, bib, notation, repro metadata). This
   changes most reported numbers, so it gates everything else.
2. Paper 1: 1.1, 1.2, 1.3 (correctness and abstract), then 1.4, 1.5
   (novelty and validation). Bring Paper 1 to submittable as the
   anchor the others cite.
3. Papers 2 to 4: complete the unfinished sections first (2.5, 3.7
   TODOs, 4.3, 4.6), then the paper-specific correctness fixes, then
   recompute every headline percentage against the Phase 0 parameters
   and correct each abstract.
4. Reposition each extension against the prior literature named above
   before reasserting novelty.

A note on scope. The four reviews did not symbolically re-derive the
r = 2 case (its Mathematica part is missing, 1.5) and reviewed the
`.Rmd` sources rather than the rendered PDFs. Items marked [S] are
suspected, not established, and should be confirmed during the fix.

---
*Rendered on 2026-06-13 at 18:22 PDT.*<br>
*Source: ~/Dropbox/prj/res/04-runin-power-analysis/runinpower/analysis/report/remediation-plan-2026-06-13.md*
