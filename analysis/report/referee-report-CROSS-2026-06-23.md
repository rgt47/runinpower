# Referee Report: Cross-Paper Consistency and Accuracy Audit
*2026-06-23 10:57 PDT*

Compendium of four manuscripts on power analysis for longitudinal
clinical trials with run-in and common close observations
(R. G. Thomas, UCSD). Reviewed as a single submitted package, with
emphasis on consistency of notation and accuracy of results *across*
the four papers, per the author's request.

- Paper 1: `01-runin-power/report.Rmd` (Woodbury closed forms)
- Paper 2: `02-correlation/report.Rmd` (AR(1) / general correlation)
- Paper 3: `03-allocation/report.Rmd` (optimal allocation)
- Paper 4: `04-gls-vs-ttest/report.Rmd` (GLS vs ANCOVA vs averaged)

---

## 1. Summary of the compendium

The four papers share a single model (a two-arm longitudinal trial
with random intercept and random slope, analysed via change scores)
and a single computational engine (the `runinpower` R package, which
forms `Sigma = R + ZGZ'` and inverts the GLS information matrix via
the Woodbury identity). Paper 1 derives closed-form treatment-effect
variances for three nested designs (no run-in; two run-in
observations; arbitrary `J_0` run-in plus `J_2` common close) under
compound symmetry and validates them against the matrix engine.
Paper 2 replaces compound symmetry with AR(1) and a five-case
correlation taxonomy. Paper 3 treats optimal allocation of a fixed
observation budget across the three phases, plus non-equal spacing,
endpoint replication, and dropout. Paper 4 compares GLS, ANCOVA, and
averaged change-score estimators under correct and misspecified
covariance. All four use a common Alzheimer's neuroimaging (MIRIAD)
numerical example.

## 2. Overall assessment and recommendation

**Major revision** (compendium-level), with the components at very
unequal maturity.

Paper 1 is close to complete and its central results are correct: I
re-derived the two-run-in closed form by hand and it is right (see
Major 1). The remaining three papers are working drafts: each
contains `<!-- TODO -->` markers, empty Discussion sections, and
results explicitly marked "to be derived." The most consequential
cross-paper problems are (i) a *displayed* derivation in Paper 1 that
does not reproduce its own (correct) boxed answer; (ii) divergent
notation and two different bibliographies across the set; (iii)
broken repository paths in Papers 2-4; and (iv) at least two
text-vs-code formula mismatches (Paper 2 cross-covariance, Paper 4
averaged estimator). None of these is fatal, but all four papers
share authoring scaffolding (`\begin{rgt}To be completed by
rgt.\end{rgt}` blocks, paired `bullets`/`orig` environments) that
must be resolved before any of them is submission-ready.

The contribution is sound and useful (a unified Woodbury treatment
that embeds the `longpower` closed forms as boundary cases and
extends them to run-in + common-close), but as a *package* the four
papers are not yet mutually consistent.

## 3. Significance and novelty

The novelty claim — no prior closed-form power formula handles
arbitrary `J_0` run-in *together with* a separate `J_2` common-close
phase under a random-slope GLS model — is plausible and survives a
literature check. The foundational reference is correctly identified:
Frost C, Kenward MG, Fox NC (2008), "Optimizing the design of
clinical trials where the outcome is a rate. Can estimating a
baseline rate in a run-in period increase efficiency?", *Stat Med*
27(19):3717-31. The positioning against `longpower`
(Iddi & Donohue, *R Journal* 2022) and the Edland/Diggle/Liu-Liang
routines is appropriate, and Paper 1's boundary-case reduction to the
Ard-Edland slope variance is the right way to establish embedding.

The synthesis (one Woodbury identity spanning four design questions)
is the real contribution and is journal-worthy for *Statistics in
Medicine*. It does not, on the evidence here, clear an *Annals*/
*Biometrika* theory bar, nor does it need to.

## 4. Major comments (correctness first)

**Major 1 — Paper 1, the displayed derivation of Eq. (var_r2k1)
does not reproduce its own boxed result, although the boxed result is
correct. [VERIFIED]**
The boxed two-run-in variance
`Var(gamma)|_{J_0=2,J_1=1} = 8 sigma^2 (sigma^2 + 5 sigma_b^2 t^2) /
[3 t^2 (sigma^2 + 2 sigma_b^2 t^2)]` is **correct** — I re-derived it
independently from `X'Sigma^{-1}X = X'R^{-1}X - W` (the boxed answer
also matches `var_gamma_r2()` in `R/variance.R`). However, the
intermediate `orig` block immediately above it
(report.Rmd ~lines 1042-1049) is wrong on two counts:
- The (3,3) cofactor is stated as `... = 6 t^4 / (16 a D)`. The
  preceding expression `1536 a D · t^4 / (16 a D)^2` actually equals
  `6 t^4 / (a D)`, not `6 t^4 / (16 a D)`. The trailing `16` is
  spurious.
- The determinant is stated as
  `det = 16 a D (3 a + 2 b t^2)^2 · t^6 / (16 a D)^3`. The true
  determinant of the combined numerator matrix is proportional to
  `a (a + 2 b t^2)(a + 5 b t^2)` (i.e. `9216 · a D (a + 2 b t^2)`
  before the `(t^2/16aD)^3` scaling). There is no `(3a + 2bt^2)^2`
  factor. Forming the cofactor/determinant ratio *as printed* yields
  `96 a D / [(3a+2bt^2)^2 t^2]`, which is **not** the boxed answer.
Why it matters: a reader who checks the algebra will conclude the
result is wrong, when it is in fact right. Remedy: replace the two
sentences with the actual cofactor `6 t^4/(aD)` and determinant
`6 t^4 (a+2bt^2)/(aD)` (so the ratio is `8aD/[3t^2(a+2bt^2)]` after
restoring the `16aD/t^2` prefactor), or simply state that the (3,3)
entry of the inverse is obtained by symbolic computation and cite
`woodbury_variance.m`.

**Major 2 — Two different bibliographies and two citation-key
conventions across the compendium. [VERIFIED]**
Paper 1 reads `01-runin-power/references.bib` (87 entries, verbose
CSL-style keys: `frostOptimizingDesignClinical2008`,
`wangTwoPeriodLinear2019`, `packerWhyHasRunIn2017`, ...). Papers 2-4
symlink `references.bib -> ../../references.bib`, i.e.
`analysis/references.bib` (68 entries, short author-year keys:
`frost2008`, `wang2019`, ...). The *same* works therefore appear
under different keys in different papers (Frost 2008 is
`frostOptimizingDesignClinical2008` in Paper 1 but `frost2008`
elsewhere). All keys do resolve within their own paper (I checked
every cited key against both files), so nothing fails to compile, but
for a compendium this is a maintenance and consistency hazard and
will produce divergent reference formatting. Remedy: adopt one shared
`.bib` and one key convention for all four papers.

**Major 3 — Broken repository paths in the Data/Code availability
and Reproducibility sections of Papers 2, 3, 4. [VERIFIED]**
Papers 2-4 state the companion manuscripts are at
`analysis/paper{1,2,3,4}-*/` and themselves at, e.g.,
`analysis/paper2-correlation/`. No `analysis/paper*` directory
exists; the actual layout is `analysis/report/0N-*/`. Paper 1 gives
the correct path (`analysis/report/{02,03,04}-*/`). Remedy: update
Papers 2-4 to the `analysis/report/` layout.

**Major 4 — Paper 2 moment-based cross-covariance: displayed formula
disagrees with the code, and contains an undefined symbol.
[VERIFIED / INSPECTED]**
The text gives
`Cov(X0_bar, XK_bar) = sigma^2/(J_0 c) · rho_1 · c^{K-1} ·
(1-c^{J_0})(1-c^{J_1})/(1-c)^2`, but the implementing chunk
`ar1_cov_phase_means()` computes
`sigma^2/(r·q) · rho · c^{K-1} · (1-c^r)(1-c^q)/(1-c)^2`. The text
divides by `J_0` only; the code divides by `J_0·J_1` (dimensionally
the code is right — it is the covariance of two phase *means*). The
symbol `rho_1` is never defined and appears to be `rho`. The same
section overloads `K` (the run-in-to-treatment gap) against `J_1`
(treatment count) and the code argument `k`. Separately,
`ar1_var_mean()` carries a factor `2*rho/(n*c)` against an `S_n`
defined with exponent `c^{i-1}`; this appears to be off by one power
of `c` relative to the standard AR(1) variance-of-mean, though I did
not fully reconcile it because the section is explicitly a draft
(`**AR(1):** (to be derived)`). Remedy: derive the section, then make
the printed formula and the code identical, and define every symbol.

**Major 5 — Paper 4 averaged-estimator definition (Eq. 2) does not
match the implemented contrast. [INSPECTED]**
Eq. (2) defines `d_ik` as the difference of a post mean over
`J_2 + 1` raw observations and a run-in mean over `J_0 + 1` raw
observations (indices `ℓ = 0..J_2` and `ℓ = 0..J_0`, i.e. baseline
included in the run-in average). The implementation `var_gamma_avg()`
applies weights `-1/r` over the `r` run-in change scores and
`1/(k+f)` over the `k+f` post change scores (baseline cancels in the
change-score representation, and is not counted among the `r` run-in
terms). The two define different estimators (the text counts
`J_0 + 1` and `J_2 + 1` observations; the code counts `r` and
`k + f`). Remedy: state the averaged estimator once, in the
change-score representation actually used by the code, and drop the
raw-`Y` indexing that does not correspond to it.

**Major 6 — Completion asymmetry across the set. [INSPECTED]**
Submitted together, the four papers are at incompatible maturity.
Paper 2 has empty Discussion, `**AR(1):** (to be derived)`, and three
`<!-- TODO -->` blocks (change-score correlations, stratified
computation, discussion). Paper 3 has empty Discussion and four TODO
blocks (threshold vs `J_1`, time-point optimisation, etc.). Paper 4
has empty Discussion, an undeveloped SUR section
(`<!-- TODO: derive the SUR variance -->`), and undelivered
"when does the gap exceed 5%" / "when to use which estimator"
sections. The simulation study in Paper 4 (5000 trials, factorial
design described in prose) is not actually run in the manuscript —
only a single-trial validation chunk executes. Abstracts promise
results the bodies do not yet deliver (e.g. Paper 2's "5-25%",
Paper 3's "5-15%"). Remedy: either complete Papers 2-4 or re-scope
their abstracts to what is delivered.

**Major 7 — Authoring scaffolding present in all four papers.
[VERIFIED]**
Every paper is built from paired `\begin{bullets}...\end{bullets}`
(a plain-language summary) and `::: {.orig}` / `\begin{orig}` (the
prose) environments, with a `\begin{rgt}To be completed by
rgt.\end{rgt}` placeholder after almost every paragraph, and Paper 1
additionally carries a standalone `\begin{rgt}To be completed by
rgt.\end{rgt}` in the Introduction. This is drafting machinery, not
manuscript content. Remedy: collapse to final prose and remove all
`rgt` placeholders before submission.

## 5. Minor comments

1. **Paper 1, Case 2 change-score listing.** "We now have `J = 3`
   change scores `(c_{i,-2}, c_{i,-1}, c_{i,1})` at times
   `(-2t, -t, 0, t)`" lists four times (including baseline `0`) for
   three change scores. Drop the `0` (it is the reference) and the
   stray `t` so the time set matches the three scores.
2. **Paper 1, Case 2 intermediate `Z`.** "`Z = t(1,1,1)' · (-2,-1,1)'
   · b_i`" is not a well-formed product (three column vectors). The
   text immediately corrects to `Z = t(-2,-1,1)'`; delete the
   malformed expression.
3. **Index conventions differ across papers.** Paper 1 uses `Y_{ij}`
   (subject, time). Papers 2 and 4 use `Y_{ijk}` (three indices,
   with group/time roles that are not stated identically). Fix one
   convention compendium-wide and define it once.
4. **"Relative efficiency" is defined with opposite orientations.**
   Paper 1's `RE = Var_base/Var_design` (`> 1` is better); Paper 4's
   `RE_avg = Var_GLS/Var_avg` (`<= 1`). Both are internally
   consistent but a shared term with inverted meaning invites reader
   error; note the orientation explicitly in each.
5. **Output formats diverge.** Paper 1 uses
   `bookdown::pdf_document2` (numbered cross-references via
   `\@ref`/`\eqref`); Papers 2-4 use plain `pdf_document`. Paper 3
   refers to "Paper 1 (Figure 1)" and "(Section 4)" as hard-coded
   text. Standardise on `pdf_document2` so cross-references resolve.
6. **`\DeclareMathOperator` sets differ** (`\Corr` absent in Paper 1,
   `\E` only in Paper 4, `\tr` declared but only used in Paper 1).
   Harmonise the preamble across papers.
7. **MATLAB filenames.** Paper 1 cites `analysis/scripts/` (correct;
   the directory holds `woodbury_variance.m`, `blup_derivation.m`,
   `raw_outcome_model.m`, `rcrm_raw_outcome.m`). The project
   `CLAUDE.md` instead names `frost311.m` / `frost312.m` under
   `analysis/data/raw_data/`; reconcile the documentation with the
   files that actually exist.
8. **Paper 3 `optimal_allocation`/`p` vs prose `J`.** Code uses `p`
   for the total budget the prose calls `J`. Align the symbol.

## 6. Missing and questionable references

| Location | Issue | Suggested action |
|---|---|---|
| All papers | Same source cited under two key schemes (`frostOptimizingDesignClinical2008` vs `frost2008`) | Unify on one `.bib` + key convention |
| Paper 3 abstract | "Frost-Kenward-Fox (2008)" author list used only here; Papers 1/2/4 say "Frost et al." | Author list is correct (verified); make attribution form consistent |
| Paper 2, AR(1) tridiagonal inverse remark | Claim that `R^{-1}` "follows from the AR(1) tridiagonal inverse" is asserted without citation; the change-score `R` is not itself tridiagonal | Add a derivation or cite a standard AR(1) precision-matrix result |
| Paper 3, Eq. `Var = 3 sigma^2/(1 - t_1 t_2) + 2 sigma_b^2` | "From the notebook derivations" with no shown derivation or unit check against Paper 1's parameterisation | Show the normalisation (`sum t_i^2 = 1`, `sum t_i = 0`) reduces to Paper 1 at the equal-spacing point |
| Paper 4 simulation | Factorial design (n, design, SNR, correlation, method; 5000 reps) described but not executed | Either run it or label it as planned |

## 7. Suggestions for strengthening

1. **Add a cross-paper consistency harness.** A single `tinytest`
   file asserting that the closed forms (Paper 1), the AR(1) limit at
   `rho -> CS` (Paper 2), the allocation objective (Paper 3), and the
   GLS variance (Paper 4) all call `var_gamma_matrix()` and agree
   numerically at shared `(J_0, J_1, J_2, sigma^2, sigma_b^2)` would
   prevent silent drift between papers.
2. **Promote Paper 1's verification chunk to a shared example.** All
   four papers re-declare the MIRIAD parameters inline
   (`0.3338 + 0.0025`, `0.9745`); expose these as a package constant
   so a single edit propagates and the four papers cannot diverge.
3. **State the embedding theorem once.** The reduction to
   Edland/Diggle/Frost (Paper 1, Eq. var_edland) is the strongest
   single result; give it a numbered proposition and reference it
   from Papers 3-4 rather than re-deriving baselines.
4. **Paper 4: a feasible-GLS (REML) simulation** rather than the
   single-trial check would substantiate the central robustness
   claim (that simpler estimators beat plug-in GLS under
   misspecification), which is currently asserted analytically via a
   sandwich formula but not demonstrated.

## 8. Scope of this review

- **Verified (re-derived or executed):** Paper 1 Case 1 variance
  `(sigma^2 + 2 sigma_b^2 t^2)/t^2`; Paper 1 Case 2 variance
  `8 sigma^2(sigma^2+5 sigma_b^2 t^2)/[3 t^2(sigma^2+2 sigma_b^2
  t^2)]` (hand re-derivation; matches `var_gamma_r2()`); the
  incorrectness of the printed cofactor/determinant intermediate
  steps; the Frost-at-`J_1=2` reduction of the Edland formula
  (`S_xx = 2t^2`); identity of MIRIAD parameters and SNR labels
  across all four papers; resolution of every cited key in each
  paper's `.bib`; the two-bibliography / two-key-scheme split;
  non-existence of the `analysis/paper*` paths cited by Papers 2-4;
  existence and contents of `analysis/scripts/` and `R/` (nine source
  files as Paper 1 states).
- **Inspected (read and judged):** the `design.R` gamma-column
  capping vs Paper 1's common-close text (consistent); `var_gamma_avg`
  / `var_gamma_ancova` vs the Paper 4 prose (the averaged-estimator
  mismatch, Major 5); the Paper 2 moment-based code vs its printed
  formula (Major 4); the completion state of each paper.
- **Not checked:** numerical reproduction of the rendered tables and
  figures (no R execution of the chunks was run); the AR(1) variance
  engine `covariance-ar1.R` and `var_gamma_ar1` against an
  independent AR(1) computation; the `optimal_allocation`,
  `common_close_beneficial`, and `var_gamma_replicated` internals;
  the sandwich-variance derivation in Paper 4's misspecification
  chunks.
- **Literature searched:** Frost/Kenward/Fox 2008 (located and
  confirmed); closed-form run-in power/sample-size formulas with
  random slopes (no prior work combining arbitrary `J_0` run-in with
  a separate `J_2` common close located, supporting the novelty
  claim); `longpower` / Iddi-Donohue 2022 (confirmed as the CRAN
  reference).

---
*Rendered on 2026-06-23 at 10:57 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/analysis/report/referee-report-CROSS-2026-06-23.md*
