# Pooled Referee Report: Run-In / Common-Close Power-Analysis Compendium
*2026-06-13 18:22 PDT*

Referee synthesis across the four manuscripts in
`analysis/report/`. Individual reports are filed alongside each
paper as `referee-report-2026-06-13.md`:

- Paper 1, `01-runin-power/` -- closed-form power under compound
  symmetry (the methodological hub)
- Paper 2, `02-correlation/` -- AR(1) and general correlation
- Paper 3, `03-allocation/` -- optimal allocation across phases
- Paper 4, `04-gls-vs-ttest/` -- GLS vs ANCOVA vs averaged
  change-score

Each paper was read in full; every headline claim was checked by
re-deriving the algebra and, where a claim was numerical, by
re-running the paper's own R chunks. The four reviews were
conducted independently and are pooled here.

---

## 1. Bottom line

The compendium rests on a sound mathematical core but is not
submittable in its present state. The shared variance engine
(Paper 1's Woodbury / Sherman-Morrison derivation and the
`var_gamma_matrix` routine) is verified correct against the Frost
(2008) closed forms and against Monte Carlo, independently, by
three of the four reviews. That is the asset worth protecting.

Around that core, however, the headline numerical claim of every
paper is false as stated, and three of the four papers are visibly
unfinished drafts with `TODO` placeholders at load-bearing
locations (including central derivations and entire Discussion
sections). Two failure modes recur across all four papers and
should be fixed once, compendium-wide, before any individual paper
is revised.

Per-paper recommendations:

| Paper | Recommendation | One-line reason |
|---|---|---|
| 1 runin-power | Major revision (borderline reject/resubmit) | Algebra sound; AD application quantitatively broken; one theorem false |
| 2 correlation | Reject in present form / major revision | Two-route equivalence unproven and self-contradictory; key formula off by a factor; headline direction backwards; unfinished |
| 3 allocation | Reject in present form / major revision | Threshold result never derived and false; novelty overstated vs optimal-design literature; unfinished |
| 4 gls-vs-ttest | Reject in present form | Two of three headline claims contradicted by own code; the Monte Carlo study does not exist; repro section fabricates seeds; unfinished |

The ordering is not accidental. Paper 1 is closest to viable
because its core is both correct and (mostly) written; Papers 2-4
inherit the correct core but add their own errors on top of an
incomplete manuscript.

---

## 2. The two compendium-wide defects (fix these first)

### 2.1 The MIRIAD / Alzheimer's parameter set is mis-specified, and the same mistake recurs in all four papers

This is the single highest-value fix. The numerical examples draw
on Frost et al.'s MIRIAD variance components but use the wrong
residual term. The papers set residual `sigma^2 = 0.0025` (Frost's
small 'additional measurement error', Table II) together with
random-slope variance `sigma_b^2 ~ 0.97`, while silently dropping
`sigma_u^2 = 0.3338`, the within-subject between-visit component
that actually plays the residual role in a change-score analysis.

Consequences observed by re-running the papers' own code:

- Paper 1: required samples collapse to n = 329 -> 2 -> 1, i.e.
  spurious 99-100% 'reductions'; Frost's own Table III gives
  175 / 122 per arm. With a sane residual (`0.3363`) Paper 1's
  reviewer recovers 385 -> 145 (a 62% reduction).
- Papers 2, 3, 4: the MIRIAD example yields N = 1-3 per group
  (Paper 4: GLS N = 1 alongside averaged N = 8193), making every
  downstream table degenerate (Paper 3's dropout strata collapse;
  Paper 4's efficiency comparison is meaningless at N = 1).

Because the signal-to-noise ratio `sigma_b/sigma` is the governing
quantity in all four papers, this single mis-specification
propagates into every headline percentage. Remedy: adopt one
correct, documented MIRIAD parameter block (residual = 0.3363 or
the explicitly justified alternative), define it once, and have all
four papers source it. This will, in turn, change most of the
reported percentages below.

### 2.2 The shared `references.bib` is broken and will silently corrupt every paper

All four papers point at one `analysis/references.bib`
(88 entries). It contains:

- `hu2021` defined two-to-three times with conflicting metadata for
  different papers sharing one key. BibTeX keeps only the last
  definition and silently drops the rest, so at least one paper
  cites the wrong source.
- Duplicate entry sets for the same work
  (`frost2008` / `frostOptimizingDesignClinical2008`;
  `wang2019` / `wangTwoPeriod...`), and inconsistent citation-key
  schemes between papers (Paper 1 uses
  `@frostOptimizingDesignClinical2008`, Paper 2 uses `@frost2008`).
- A corrupted Laird-Ware entry; a possibly wrong `winkens2005` DOI;
  roughly 50 uncited Zotero-dump entries.

Remedy: de-duplicate to a single canonical key per work, fix the
`hu2021` collision first (it is an active mis-citation, not
cosmetic), prune uncited entries, and standardise keys across the
compendium.

### 2.3 Secondary shared issues

- Notation drift `J0/J1/J2` (prose) vs `r/k/f` (code and R
  package) runs through all four papers. Symbol collisions
  compound it: `a,b` are used both for subject random effects and
  for `sigma^2, sigma_b^2`; Paper 2 redefines `c = rho` whereas
  Paper 1 and the Mathematica sources use
  `c = sigma^2 - sigma_b^2 t^2`. Choose one notation and a symbol
  table, apply across the set.
- Reproducibility paths are wrong project-wide: papers cite
  `analysis/paperN-*/`; the actual layout is
  `analysis/report/0N-*/`.
- Derivation-source drift: CLAUDE.md and the manuscripts reference
  `frost311.m` / `frost312.m`, but the actual Mathematica files are
  `woodbury_variance.m`, `blup_derivation.m`, `raw_outcome_model.m`,
  `rcrm_raw_outcome.m`. The scripts are Mathematica/Wolfram but are
  called 'MATLAB' in the text. Note also that `woodbury_variance.m`
  Part 2 derives the `r = 1` case, not Paper 1's headline `r = 2`
  case, so the symbolic trail for `var_r2k1` is currently missing.

---

## 3. What is verified correct (the protected core)

Stated plainly because it matters for triage: the foundation is
sound, so the path forward is correction and completion, not
re-derivation.

- Paper 1's `var_r0k2` matches Frost Eq. 13; `var_gamma_matrix(1,1,0)`
  matches Frost Eq. 16 exactly across 12 parameter sets; `var_r2k1`
  matches the matrix routine to 3.6e-15; the Edland recovery
  identity is exact; Sherman-Morrison and Woodbury are correctly
  applied, including the per-subject vs summed outer-product
  subtlety. [VERIFIED, Paper 1 review]
- The same `var_gamma_matrix` engine and the closed form
  `(sigma^2 + 2 sigma_b^2 t^2)/t^2` reappear unchanged and correct
  in Paper 3 and Paper 4, each cross-checked against Monte Carlo.
  [VERIFIED, Papers 3 and 4 reviews]
- Paper 2's GLS routine `var_gamma_ar1` matches an independent
  block-diagonal GLS to machine precision across `J0` and `rho`;
  the AR(1) change-score covariance formula is also correct.
  [VERIFIED, Paper 2 review]
- Paper 4's estimator-ordering claim GLS <= ANCOVA <= averaged
  holds with zero violations across 413 grid points (the inequality
  is true even though the paper only asserts, never proves it).
  [VERIFIED, Paper 4 review]

The cross-paper agreement on the core formula is itself evidence of
good internal discipline. The problems are concentrated in the
application layer, the prose claims, and the unfinished sections.

---

## 4. Per-paper headline-claim failures

All entries below were confirmed by re-running each paper's own R
chunks unless marked otherwise.

### Paper 1 (runin-power)

1. AD application misuses the MIRIAD components (see 2.1); sample
   sizes are nonsensical (329 -> 2 -> 1). [VERIFIED]
2. The claim that `J0 = 1` worsens power relative to `J0 = 0` is
   FALSE. An exhaustive grid shows `J0 = 1` is never worse, and
   exactly equal at `sigma_b^2 = 0.5, t = 1, k = 2`. Correct
   statement: 'never worse, sometimes no gain.' The Frison-Pocock
   citation does not support the original claim. [VERIFIED]
3. Abstract '30-45% reduction' is unsupported; it is Frost's MIRIAD
   headline, not this paper's result (the paper's own run-in table
   gives 14 / 28 / 76%). [VERIFIED]

### Paper 2 (correlation)

1. The two-route (Woodbury vs moment-based) equivalence, the
   paper's headline internal-validation argument, is NOT
   demonstrated and is self-contradictory: the abstract claims
   identical results, the only numerical comparison prints
   disagreeing values (1.4375 vs 1.6662), and Section 5 states the
   routes 'compute different quantities.' They target different
   estimands. [VERIFIED]
2. The moment-based variance-of-mean formula is wrong: the lag sum
   uses `c^{i-1}` where `c^{i}` is required, off by a factor `c`
   per term (n = 2, rho = 0.3: 1.000 vs correct 0.650). [VERIFIED]
3. The central '5-25% higher sample size under AR(1)' claim is
   false in both direction and magnitude: AR(1) sizes are lower
   than or equal to compound symmetry, and the deviation grows with
   rho, so 'largest at small rho' is also backwards. [VERIFIED]
4. Novelty premise is false: the paper claims `longpower` has 'no
   analytic AR(1)-residual formula', but `power.mmrm.ar1` is
   exactly that (Lu-Luo-Chen 2008). [VERIFIED]

### Paper 3 (allocation)

1. The common-close benefit threshold, the abstract's lead result,
   does not exist: `common_close_beneficial()` returns TRUE for all
   `k` and all ratios down to 0.001; the reported 'threshold 0.1'
   is merely the grid floor. A `TODO` admits it was never derived.
   [VERIFIED]
2. Run-in is never optimal under fixed budget: optimal `J0* = 0` in
   every cell of a wide sweep, contradicting the paper's framing
   and intended message. [VERIFIED]
3. The non-equal-spacing variance formula is model-inconsistent: it
   is derived from a raw-outcome, random-slope-only model
   (`raw_outcome_model.m`) while the rest of the paper uses change
   scores with `R = sigma^2(I + 11')`; it does not match direct
   GLS. [VERIFIED mismatch; SUSPECTED error]
4. The 'optimisation' is a brute-force grid search presented as a
   'closed-form decision procedure'; no KKT or monotonicity
   analysis is given. [VERIFIED]

### Paper 4 (gls-vs-ttest)

1. Abstract's central 'ANCOVA efficiency loss below 5% across a
   broad regime' is false: re-running the paper's grid,
   RE(GLS/ANCOVA) is mostly 0.1-0.5 (ANCOVA variance 2-10x GLS);
   the >= 0.95 condition holds in only 5.2% of cells, all with
   `J1 = 1`. The formula is correctly implemented; the verbal claim
   is wrong. [VERIFIED]
2. The robustness claim ('simpler estimators can outperform GLS
   under misspecification') is contradicted by the paper's own
   figures: GLS retains a roughly 10x variance advantage even when
   misspecified. [VERIFIED]
3. The Monte Carlo study does not exist. The text claims a
   5000-replicate factorial simulation with empirical variance,
   power, and coverage; the `simulation` chunk defines functions
   and runs ONE trial under `set.seed(42)`. There is no replicate
   loop, no factorial sweep, no power or coverage, no MC standard
   errors. The Reproducibility section additionally cites a seed
   `set.seed(20260418)` with per-replicate seeding that appears
   nowhere in the source, and names `nlme` as the GLS fitter though
   `nlme` is never called. This is a fabricated-methods problem,
   not a coding bug, and must be corrected before resubmission to
   any journal. [VERIFIED]

---

## 5. Novelty assessment of the set

Pooled verdict: the compendium's genuine, defensible contribution
is narrower than claimed and is concentrated in Paper 1.

- The defensible novel element across the set is the closed-form
  treatment of the common-close phase under the random-slopes
  model, and its integration with run-in averaging in a single
  Woodbury-based framework. That is real and, in Paper 1, mostly
  correct, though it is currently delivered without a closed form
  for the `J2 > 0` case and without numerical validation (no test
  exercises `J2 > 0`).
- The `J0` generalisation is Woodbury bookkeeping on Frost (2008),
  which already supplies the general LMM machinery, the run-in
  design, and the MIRIAD application.
- Paper 2's AR(1) extension is undercut by existing analytic AR(1)
  tooling (`longpower::power.mmrm.ar1`) and by the
  summary-statistics literature (Frison-Pocock 1992, in the bib but
  uncited; Munoz et al. 1992, missing; Winkens et al. 2007).
- Paper 3 overstates novelty against a well-developed optimal-design
  literature it does not cite: Tekle-Tan-Berger, Fedorov-Hackl /
  Fedorov-Jones, ODmixed (optimal design under dropout), and a 2023
  BMC Medical Research Methodology paper on optimal pre-post
  allocation that overlaps its run-in-vs-replication material. Only
  Winkens (2005/2007) from this stream is cited.
- Paper 4 reproduces the Frison-Pocock (1992/1997) POST / CHANGE /
  ANCOVA efficiency ordering; those papers are in the shared bib but
  are not cited in Paper 4, the paper they bear on most directly.

Taken together, the set currently reads as one solid methodological
note (Paper 1, after repair) plus three extensions that are either
incremental or not yet positioned against their prior literature.

---

## 6. Recommended remediation order

1. Fix the shared MIRIAD parameter block (2.1) and the shared
   `references.bib` (2.2). These two touch every paper and change
   most reported numbers.
2. Standardise notation, citation keys, repro paths, and the
   derivation-source naming across the compendium (2.3).
3. Repair Paper 1's two false/overstated claims (the `J0 = 1`
   theorem and the abstract percentage), add a closed form and a
   validation test for the `J2 > 0` common-close case, then bring
   Paper 1 to a submittable state as the anchor.
4. For Papers 2-4: complete the unfinished sections, then correct
   the paper-specific errors (Paper 2 moment formula and
   two-route argument; Paper 3 threshold derivation and
   spacing-model consistency; Paper 4 the actual Monte Carlo study
   and the Reproducibility section). Re-derive each headline
   percentage after step 1, since several will change sign or
   magnitude.
5. Position each extension against the prior literature named in
   Section 5 before reasserting novelty.

---

## 7. Scope of this review

- Verified (re-derived or executed): the core variance algebra in
  all four papers; every headline numerical claim, by re-running
  the papers' own R chunks; the Paper 2 moment-formula error and
  the absence of the Paper 4 Monte Carlo loop; the MIRIAD
  parameter consequences.
- Inspected (read and judged): manuscript structure, abstract
  accuracy, completeness (TODO placeholders), citation-to-claim
  support, bibliography hygiene, and the Mathematica derivation
  sources for the available cases.
- Not checked: full symbolic re-derivation of the `r = 2` case
  (the corresponding Mathematica part is missing); exhaustive
  re-execution of every chunk in every paper; the rendered PDFs
  (only the `.Rmd` sources were reviewed).
- Literature searches run: analytic AR(1) power tooling and serial
  -correlation longitudinal methods; optimal-design literature for
  repeated-measures and mixed models; the ANCOVA-vs-change-score
  efficiency literature; run-in and delayed-start design methods;
  Frost / Kenward / Fox and the MIRIAD application. Details in each
  paper's individual report, Section 8.

---
*Rendered on 2026-06-13 at 18:22 PDT.*<br>
*Source: ~/Dropbox/prj/res/04-runin-power-analysis/runinpower/analysis/report/referee-report-POOLED-2026-06-13.md*
