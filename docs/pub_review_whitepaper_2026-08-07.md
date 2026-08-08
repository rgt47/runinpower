# Publication Review White Paper: Run-In Power Analysis Compendium
*Review date: 2026-08-07 18:05 PDT*

This white paper reports a full journal-referee-style re-review of the
four manuscripts in `analysis/report/`, conducted against the standards
of Statistics in Medicine, Biometrics, and comparable venues. Each
manuscript was read in full by an independent review pass, headline
numerical claims were re-executed against the `runinpower` package
source (`R/*.R`), and every finding from the 2026-06-13 referee reports
and the 2026-06-23 cross-paper audit was re-checked for resolution
status. Epistemic labels used throughout: verified (code executed),
inspected (source read), inferred, unverified.

A repository action taken during this review session: all 119
`\begin{rgt}To be completed by rgt.\end{rgt}` placeholder blocks were
removed from the four `report.Rmd` files at the author's request
(48 / 17 / 31 / 23 in Papers 1-4 respectively; verified zero remain).
The findings below describe the manuscripts as reviewed; the
bullets/`.orig` two-layer drafting scaffold remains in place and is
still a submission blocker (Major Issue 1).

---

## 1. Summary of the work under review

**Paper 1, `01-runin-power/report.Rmd` (2,039 lines).** Derives the
variance of the GLS treatment-effect (slope-difference) estimator in a
two-arm longitudinal trial under a random-slopes linear mixed model
with change-score parameterization, via the Woodbury identity. Three
nested cases: no run-in (recovering Frost et al. 2008 and, as boundary
cases, the Ard-Edland and Diggle formulas underlying `longpower`); a
closed form for two run-in observations; and a general
(J0, J1, J2) configuration including a post-treatment common-close
phase evaluated by the package's matrix routine. The MIRIAD
Alzheimer's example now uses the corrected residual
(sigma^2 = 0.3338 + 0.0025 = 0.3363) and yields per-group sample
sizes 385 (standard) to 145 (two run-in, 62%) to 69 (full design,
82%), all reproduced from code (verified). The quantitative core is
sound; the manuscript around it is not finished.

**Paper 2, `02-correlation/report.Rmd` (829 lines).** Replaces
compound symmetry (CS) with AR(1) residual correlation, presenting a
Woodbury/GLS route (`var_gamma_ar1`, verified correct against
brute-force GLS) and a moment-based phase-mean route. Since June the
moment-based variance and cross-covariance formulas have been fixed
(verified), and the MIRIAD units repaired. However, the abstract's
headline claim, that AR(1) inflates required sample sizes by 5-25%
with the largest deviation at small rho, is contradicted by the
paper's own functions: AR(1) sample sizes are 4-36% *smaller* than CS
at rho = 0.5, and the dominant deviation grows with large rho
(verified). Two core derivations still read "(to be derived)."

**Paper 3, `03-allocation/report.Rmd` (1,268 lines).** Poses a
discrete allocation problem: split a fixed budget of J post-baseline
observations across run-in, treatment, and common close to minimize
the treatment-effect variance, with extensions to non-equal spacing,
endpoint replication, and dropout. The verified central computation
shows the optimum places *zero* run-in observations in every cell of
the parameter sweep, yet the newly drafted Discussion asserts that
"one or two run-in observations" are optimal, an outright internal
contradiction with the paper's own Table 2 (verified). The abstract's
threshold claim for common-close benefit is a grid-floor artifact
(verified: benefit holds down to ratio 0.001; no threshold exists).

**Paper 4, `04-gls-vs-ttest/report.Rmd` (1,309 lines).** Compares GLS
(known covariance), ANCOVA on phase means, and an averaged
change-score t-test estimator, with an analytic sandwich
misspecification analysis and a Monte Carlo study (new since June:
1,000 replicates, seeded, real). The verified ordering
GLS <= ANCOVA <= averaged holds without exception across the paper's
grid, but both abstract headlines remain false against the paper's own
code: the "ANCOVA loss below 5%" claim (actual RE is 0.29-0.36 at the
representative designs; the <5% condition holds only in degenerate
J1 = 1 cells) and the "simpler estimators can outperform GLS under
misspecification" claim (GLS retains a large margin in every sandwich
sweep and every simulation cell; verified).

**Coherence across the set.** The four papers share a verified-correct
computational engine and a common MIRIAD example, which is a genuine
strength. They diverge in bibliography keys (Paper 1 uses its own
87-entry `references.bib` with CSL-style keys; Papers 2-4 symlink the
shared 66-entry file with short keys), in output format
(`bookdown::pdf_document2` vs plain `pdf_document`), in notation
(J0/J1/J2 in prose vs r/k/f in all code, a/b symbol collisions, Paper
2's c = rho vs Paper 1's c = sigma^2 - sigma_b^2 t^2), and in the
orientation of "relative efficiency." Papers 2-4 cross-reference
"Paper 1" as internal compendium text that cannot survive standalone
submission. There is no material redundancy problem; the division of
material across the four papers is clean. The coherence problem is
consistency of claims: Paper 1 advertises run-in gains while Paper 3's
own optimizer shows zero run-in is optimal under a fixed budget, a
tension no paper acknowledges (the reconciliation, fixed-budget vs
added-visit framings, is exactly what Paper 3 should be about).

## 2. Major issues

1. **All four manuscripts remain drafting scaffolds, not papers.**
   (All papers, global.) Every paragraph exists as a `bullets` skeleton
   plus an `.orig` prose block (the 119 `rgt` placeholders were removed
   during this session, but the two remaining layers still render as
   annotated drafting apparatus, with colored environments and
   paragraph counters via `claudecode.tex`). No journal will referee a
   document in this form. Remediation: collapse each paper to
   single-layer final prose and delete the scaffold environments and
   preamble machinery.

2. **Paper 1 retains a false theorem previously disproved.**
   (`01-runin-power/report.Rmd`, Section 4, approx. lines 1515-1540.)
   The claim that one run-in observation can be worse than none, with
   the design recommendation "plan either zero or at least two, never
   one," was demonstrated false in the 2026-06-13 report and was
   re-verified false in this review (0 of 750 configurations worse;
   verified). A false claim surviving a documented disproof is the most
   serious category of finding. Remediation: replace with the
   never-worse, sometimes-no-gain characterization, derive the
   zero-gain boundary, and drop or re-scope the Frison-Pocock citation.

3. **Paper 2's two headline claims are contradicted by its own code.**
   (`02-correlation/report.Rmd`, abstract and Discussion.) (a) The
   5-25% AR(1) inflation claim is wrong in direction, magnitude, and
   localization (verified: AR(1) sample sizes are smaller at matched
   nominal rho except in a narrow J0 = 4, low-rho corner of about 5%).
   (b) The "both routes yield identical results, providing internal
   validation" claim is contradicted by Section 5 and by the chunk
   output (1.2188 vs 1.6662; verified); the routes target different
   estimands. Remediation: recompute and re-sign every quantitative
   claim; state the cross-structure matching convention; withdraw or
   prove the equivalence claim (which first requires fixing
   `var_gamma_ar1` at sigma_b^2 = 0, where `solve(solve(G))` is
   singular; verified failing).

4. **Paper 3's Discussion contradicts its own Results.**
   (`03-allocation/report.Rmd`, chunk `optimal-grid` vs Discussion.)
   The optimizer returns J0* = 0 and J2* up to 4 in every cell
   (verified), while the Discussion asserts 1-2 run-in and at most one
   common close, and the Introduction claims the first run-in
   observation yields the largest gain, which Section 3.2 itself
   denies. Additionally the abstract's threshold claim is a grid-floor
   artifact (verified), and "closed-form decision procedure" describes
   an exhaustive grid search (inspected). Remediation: make the
   zero-run-in-under-fixed-budget result the honest headline,
   reconcile with Paper 1's added-visit framing, and rewrite abstract,
   Introduction, and Discussion to agree with the computed tables.

5. **Paper 4's abstract is refuted by every exhibit meant to support
   it.** (`04-gls-vs-ttest/report.Rmd`, abstract, Sections 4-5.) The
   "below 5% ANCOVA loss" claim fails (verified RE 0.29-0.36 at
   designs of interest; the MIRIAD table shows N_ancova up to roughly
   13 times N_GLS), the promised "conditions characterizing this
   regime" are an empty TODO, and the robustness claim is contradicted
   by both sandwich sweeps and the paper's own simulation table
   (verified). Remediation: reframe around what is verified (GLS
   dominates; the simplicity penalty is large and characterizable in
   closed form), or introduce the estimator for which a near-GLS claim
   could hold (per-visit ANCOVA or cLDA) and re-derive.

6. **Validation and reproducibility statements claim work that does
   not exist.** (Papers 1, 2, 3; Paper 4 partially repaired.) Papers
   1-3 cite a `set.seed(20260418)` Monte Carlo benchmark that appears
   nowhere (inspected); Paper 4's simulation now exists but its
   Reproducibility section describes per-replicate `+ rep_idx` seeding
   not present in the code, names `nlme` as the GLS fitter though
   `nlme` is never called, and delivers a fraction of its stated
   factorial design with no coverage or size control and no Monte
   Carlo standard errors (verified). These are fabricated-methods
   findings and are mandatory fixes before any submission.
   Remediation: implement the described validation (cheap: the Paper 1
   reviewer's direct-GLS and Monte Carlo checks each ran in seconds
   and passed) or delete every claim to it.

7. **The novel contribution of Paper 1 is still asserted, not
   derived.** (`01-runin-power/report.Rmd`, common-close section.) The
   frozen-offset construction for the common-close phase presumes zero
   carryover and perfect level maintenance, precisely where
   Alzheimer's disease-modification debates live; no estimand
   statement, no assumptions, no contrast with delayed-start designs
   (Liu-Seifert 2015 is in the bibliography, uncited). The r = 2
   closed form and the general X'R^{-1}X entries are numerically
   correct (verified to machine precision) but have no symbolic
   derivation trail; `woodbury_variance.m` covers r = 1 only
   (verified). Remediation: derive the common-close mean from the
   generative model, state assumptions, add a Mathematica script or
   appendix for the r = 2 and general cases.

8. **Literature positioning is materially incomplete or incorrect in
   Papers 2-4.** Paper 2 claims an analytic AR(1) formula gap that
   `longpower::power.mmrm.ar1` (Lu-Luo-Chen 2008) already fills
   (`lu2008` sits in the bibliography unused; inspected), and
   miscasts the rho -> 0 limit as a "compound-symmetry limit" of
   AR(1). Paper 3 cites none of the optimal-design literature
   (Tekle-Tan-Berger, Fedorov, Atkinson-Donev, Mentre; verified absent
   from the bibliography). Paper 4 never cites Frison and Pocock
   (1992, 1997), the closest prior work on multi-baseline
   ANCOVA/change-score orderings, though both entries are in the
   bibliography, and omits Liang-Zeger 1986 at the sandwich
   construction (inspected). Remediation: engage each and restate the
   narrower defensible novelty claims.

## 3. Minor issues

1. `zhao2022` is duplicated verbatim in the shared
   `analysis/references.bib` (lines 505 and 758; verified), a
   regression introduced since the June cleanup; Paper 1's separate
   bibliography and key scheme still diverge from Papers 2-4.
2. British spellings persist in abstracts and prose of all four papers
   (randomisation, characterised, minimise), inconsistent within
   documents and contrary to the project's US-English standard.
3. Notation: J0/J1/J2 prose vs r/k/f code with no stated mapping in
   any paper; a/b collide with the random effects a_i, b_i; Paper 2's
   c = rho vs Paper 1's c reparameterization.
4. Paper 1: three-change-score time list includes four times; the
   malformed Z product expression; circular H display; literal `J_0`
   axis labels in `filled.contour`; the AD table's filter hides the
   baseline row that the 62%/82% reductions are computed against;
   "MATLAB" for what are Mathematica scripts; hard-coded
   `/Users/zenn/Dropbox` path in the footer chunk; `nash2021`
   (slopepower), `galbraith2002`, `liuseifert2015` in bibliography but
   uncited.
5. Paper 2: robustness figure `ylim = c(0.5, 2)` clips its own curves
   (values reach 0.184; verified); moment-based functions live only in
   the Rmd, contradicting the Data Availability claim; `geometry
   right=5cm` drafting artifact; K convention under-specified.
6. Paper 3: run-in-vs-replication dominance claim fails at ratio = 1
   with one extra measurement (verified); `var_gamma_replicated`
   contains dead code and a vacuous branch (inspected); raw `cat()`
   output presented as results; `design.R` misdescribed as an
   allocation routine.
7. Paper 4: "MVUE" should be BLUE; z rather than t critical values
   despite the t-test framing; "unaffected" captions contradicted by
   the plotted rho-dependence; leftover `set.seed(42)` debug output;
   fragile cross-chunk variable coupling (`tt`, silently redefined
   `a`); title omits ANCOVA, the estimator the paper actually
   recommends.
8. Papers 2-4 cite "tinytest" in Reproducibility while the Makefile
   and CI use `devtools::test()`/testthat; the repository in fact
   contains both `tests/tinytest.R` and `inst/tinytest/` (inspected),
   so the Makefile documentation, not the manuscripts, may be the
   stale side; reconcile once, project-wide.
9. Abstracts in Papers 1-4 reference "Paper 1 ... Paper 4 in this
   compendium"; replace with citable preprint references at
   submission. No paper gives a repository URL or DOI in its Data
   Availability statement.

## 4. What remains to be done

Ordered by importance for submission readiness.

**(a) Required for correctness**

1. Paper 1: replace the false J0 = 1 theorem (Major 2).
2. Paper 2: recompute and re-sign the AR(1)-vs-CS comparison; withdraw
   or prove the two-route equivalence; fix `var_gamma_ar1` at
   sigma_b^2 = 0 (Major 3).
3. Paper 3: resolve the Results-vs-Discussion contradiction; replace
   the threshold artifact with the genuine fixed-budget trade-off;
   drop "closed-form" for the grid search (Major 4).
4. Paper 4: delete or repair both abstract headlines; either prove the
   ordering proposition (two lines: both competitors are linear
   unbiased, GLS is BLUE; ANCOVA optimizes over phi where averaged
   fixes phi = 1) or stop asserting the boundary (Major 5).
5. All papers: remove or implement every claimed-but-absent validation
   (phantom Monte Carlo seeds, `nlme`, per-replicate seeding)
   (Major 6).
6. Paper 3: rebuild the dropout strata as genuine monotone-dropout
   patterns; reconcile or cut the non-equal-spacing section, which
   lives in a different model (Section 6.1).

**(b) Required for acceptance**

7. Collapse the two-layer scaffold to final prose in all four papers
   and delete the `claudecode.tex` drafting machinery (Major 1).
8. Paper 1: derive the common-close estimand with stated carryover
   assumptions and a delayed-start contrast; supply the r = 2 and
   general-entry symbolic derivations (Major 7).
9. Paper 4: run the full stated simulation (factorial design, null
   case for size, coverage, MC standard errors, AR(1)-truth arm,
   feasible-GLS and feasible-ANCOVA arms; cache results to RDS).
10. Literature positioning per paper (Major 8): Paper 2 vs
    `power.mmrm.ar1` and Frison-Pocock 1992; Paper 3 vs the
    optimal-design literature; Paper 4 vs Frison-Pocock and
    Liang-Zeger; add cLDA to Paper 4 or justify its exclusion.
11. Unify the bibliography (one file, one key scheme, de-duplicate
    `zhao2022`, prune uncited entries) and the notation (one symbol
    table; state the J/r-k-f mapping; resolve the a/b and c
    collisions; single RE orientation).
12. Complete or cut every remaining "(to be derived)"/TODO section
    (Paper 2 change-score correlations and stratified computation;
    Paper 3 threshold and spacing sections; Paper 4 SUR, 5% boundary,
    decision-rule sections).

**(c) Desirable polish**

13. Standardize US English; standardize on `bookdown::pdf_document2`
    so cross-references resolve; remove `geometry right=5cm`.
14. Add a cross-paper consistency tinytest asserting the four papers'
    engines agree at shared parameter points; expose the MIRIAD
    parameter block as a package constant sourced by all four papers.
15. Replace compendium-internal "Paper N" references with citations;
    add repository URL/DOI to Data Availability; small-sample
    (Kenward-Roger) caveat where n is small; figure-quality passes
    (plotmath labels, unclipped axes, no raw `cat()` output).

## 5. Recommended framing

**Paper 1.** Plausible framings: (i) general closed-form power
formulas for run-in designs (current); (ii) the common-close phase as
the headline, with run-in as the established base; (iii) a
software-plus-embedding paper unifying the `longpower`/`slopepower`
closed forms. Recommendation: (ii). The J0 generalization is Woodbury
bookkeeping on Frost (2008) and overlaps the author's own Hu-Mackey-
Thomas (2021) line, whereas no located prior work combines arbitrary
run-in with a separate common-close phase under a random-slope GLS
model (literature check in the 2026-06-23 audit and re-affirmed this
round; inferred, not exhaustive). Implications: title and abstract
lead with the common close; comparators are Frost 2008, Wang 2019,
Liu-Seifert 2015 (delayed start), `longpower`, `slopepower`; target
Statistics in Medicine or Pharmaceutical Statistics. Emphasize the
embedding proposition (reduction to Edland/Diggle/Frost) as a numbered
result; move the r = 2 cofactor algebra to an appendix or script
citation.

**Paper 2.** Plausible framings: (i) AR(1) inflation with a five-case
taxonomy (current); (ii) robustness audit: how wrong is the CS-based
run-in formula when residual correlation is AR(1); (iii) GLS vs
phase-mean estimator efficiency under AR(1). Recommendation: (ii).
Framing (i) is untenable because the inflation claim is false under
the paper's own model and the taxonomy is empty; framing (iii)
duplicates Paper 4's theme. The honest, verified result is that CS is
generally conservative at matched nominal rho with a narrow
anticonservative corner, plus the one solid novel object: the
closed-form AR(1) GLS variance for run-in/common-close random-slope
designs, which `power.mmrm.ar1` does not cover. Implications: retitle
toward robustness of run-in power formulas to serial correlation;
comparators are Lu-Luo-Chen 2008, Frison-Pocock 1992, Winkens 2007;
target Pharmaceutical Statistics or Statistical Methods in Medical
Research; the five-case taxonomy moves to a short subsection or is
dropped.

**Paper 3.** Plausible framings: (i) optimal allocation with
closed-form thresholds (current); (ii) zero run-in under a fixed visit
budget, reconciling run-in efficiency with allocation optimality;
(iii) numerical design guidance for three-phase AD trials.
Recommendation: (ii). The counterintuitive verified finding, that
run-in visits are never chosen when visits are fungible while Paper
1's gains arise only when run-in visits are additional, is the
genuinely interesting result already sitting unclaimed in the paper's
own output, and it converges with the pre-post allocation literature
rather than competing with it. Implications: retitle around the
reconciliation; the common-close allocation (J2* growing with
sigma_b/sigma) becomes the secondary result and the defensible version
of the threshold story; add the optimal-design citations; target
Statistics in Medicine or Pharmaceutical Statistics; the replication
and spacing extensions move to supplementary material unless
completed.

**Paper 4.** Plausible framings: (i) ANCOVA is nearly as good as GLS
and more robust (current); (ii) quantifying the efficiency cost of
simple analyses (change-score, ANCOVA-on-means) in run-in trials;
(iii) when do simple estimators beat feasible GLS. Recommendation:
(ii), with one cell of (iii) added. Framing (i) is refuted by the
paper's own verified numbers on both halves. Everything framing (ii)
needs already exists and verifies: the strict ordering, closed-form
relative efficiencies in the run-in parameterization, penalties of
3-10x growing with J1 and signal-to-noise ratio, and MIRIAD
consequences. The one genuine robustness question (small-sample
feasible GLS with REML-estimated components, where crossovers
plausibly do occur) deserves a dedicated simulation section rather
than the current false claim. Implications: retitle to name all three
estimators and the penalty question; position against Frison-Pocock;
drop SUR; target Statistics in Medicine or Pharmaceutical Statistics
as a compact methods paper.

## 6. Assessment

Compendium-level verdict: **major revision, not currently evaluable
for acceptance in any component.** Per-paper: Paper 1, major revision
(quantitative core sound and verified, but a disproved theorem
persists, the novel contribution is underived, and the validation
claims are phantom); Papers 2 and 3, reject in present form with
encouragement to resubmit (trustworthy computational engines beneath
abstracts and discussions contradicted by their own code); Paper 4,
reject in present form pending reframing (the honest paper the
verified numbers support is the opposite of the one the abstract
claims).

The asset worth protecting is real and has been independently verified
three times now: the Woodbury variance engine, the corrected MIRIAD
example, and the strict estimator ordering are correct to machine
precision against direct GLS and Monte Carlo. The compendium's failure
mode is equally consistent: headline claims written before, and never
reconciled with, the computations beneath them. The June-to-August
trajectory shows genuine progress on parameters, bibliography hygiene,
paths, and Paper 4's simulation, but near-zero progress on the
claims-level corrections that all three prior reports identified as
dispositive, and two papers now contain internal contradictions that
any referee will find on a first read. The fastest path to a
submission is the one the remediation plan already prescribes: fix the
claims to match the verified code, not the reverse, starting with
Paper 1 as the anchor.

## 7. Revision history

- 2026-08-07: Initial pub_review white paper. Baseline drawn from the
  2026-06-13 per-paper referee reports, the pooled report, the
  remediation plan, and the 2026-06-23 cross-paper audit and numerical
  verification. Resolved since June: MIRIAD parameter fix in all four
  papers (verified), `hu2021` bibliography de-duplication, repository
  path corrections, Paper 2 moment-based formula fixes, Paper 4
  simulation now exists, Discussions drafted in Papers 2-4. Newly
  identified: Paper 3 Discussion-vs-Results contradiction; Paper 2
  robustness figure clipping; Paper 4 power-column anticonservatism
  artifact and inaccurate seeding description; `zhao2022` duplicate
  bibliography key; two-layer drafting scaffold across all papers.
  Still open: Paper 1 false J0 = 1 theorem, common-close derivation,
  r = 2 symbolic trail; Paper 2 headline direction and equivalence
  claims; Paper 3 threshold artifact and closed-form overclaim; Paper
  4 both abstract headlines; phantom validation claims in Papers 1-3.
  Repository action this session: removed all 119 "To be completed by
  rgt." placeholder blocks from the four manuscripts.
- 2026-08-08: Remediation pass (not a fresh review). Applied and
  verified against the package source: Paper 1 false J0 = 1 theorem
  replaced, phantom Monte Carlo replaced by a real seeded validation
  chunk (direct GLS vs Woodbury, 396 configurations including f > 0),
  abstract superlative removed, MATLAB/Mathematica and mechanical
  fixes. Paper 2 abstract and Discussion re-signed to the verified
  direction (AR(1) generally smaller than CS at matched nominal rho;
  anticonservative corner at most about 6%), two-route equivalence
  claim withdrawn and non-coincidence at sigma_b^2 = 0 documented,
  CS-limit and longpower positioning corrected, Section 2.2
  change-score convention fixed and the AR(1) correlation derivation
  completed and numerically verified. Paper 3 rewritten around the
  verified zero-run-in headline, threshold artifact replaced by a
  computed fixed-budget crossover table, closed-form language removed,
  J0 = 1 anomaly restricted to its verified band (0.53-1.51), dropout
  strata rebuilt as genuine monotone-truncation patterns. Paper 4
  abstract re-signed (GLS dominates; RE >= 0.95 essentially only at
  J1 = 1), ordering proposition stated and proved, gap
  characterization chunk added, misspecification captions and
  simulation design description made honest, seeding/nlme
  misstatements fixed. Shared fixes: var_gamma_ar1 rewritten in
  push-through form (valid at sigma_b^2 = 0; all 61 package tests
  pass), duplicate zhao2022 bibliography entry removed. Papers 2-4
  additionally pass an end-to-end purled-chunk execution. Newly
  identified during remediation (open): the averaged and ANCOVA
  contrasts in Paper 4 estimate c times gamma with
  c = t(J1(J1+1)/2 + J1 J2)/(J1+J2), so the raw tabulated penalties
  overstate the equal-estimand penalty by c^2 (up to about 5x at
  (4,3,1)); disclosed in the manuscript, rescaled tables deferred to
  the author. Still open (author-level): common-close estimand
  derivation and r = 2 symbolic trail (Paper 1), taxonomy content and
  stratified-variance TODO (Paper 2), Section 6 spacing model and
  optimal-design citations (Paper 3), cLDA/SUR/feasible-GLS
  simulation and Frison-Pocock citations (Paper 4), scaffold collapse
  to final prose and compendium-internal cross-references (all).

---
*Rendered on 2026-08-07 at 18:05 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/docs/pub_review_whitepaper_2026-08-07.md*
