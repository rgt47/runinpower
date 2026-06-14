# Referee Report: 'GLS versus Averaged Change-Score Estimators for Clinical Trials with Run-In Observations'
*2026-06-13*

Manuscript: Paper 4 of the run-in power-analysis compendium (Thomas, UCSD).
Reviewed against: the .Rmd source, the package source it calls
(`R/variance.R`, `R/ancova.R`, `R/averaged.R`, `R/design.R`,
`R/covariance.R`, `R/covariance-ar1.R`), the shared bibliography
(`analysis/references.bib`, symlinked as `references.bib`), and the
Mathematica derivations in `analysis/scripts/`. All numerical claims were
re-executed against the package code (R 4.5.3, `devtools::load_all`).

---

## 1. Summary

The paper compares three estimators of a longitudinal treatment effect in
trials carrying run-in observations: (i) the GLS estimator under a
random-slopes mixed model, taken from Paper 1; (ii) an ANCOVA estimator
that regresses the post-randomisation mean on the run-in mean; and (iii) an
averaged change-score estimator analysed by a two-sample t-test. The stated
contributions are closed-form relative-efficiency (RE) expressions for the
three estimators as functions of the design $(J_0, J_1, J_2)$ and the
variance ratio $\sigma_b/\sigma$, an analytic robustness assessment under
misspecified $G$ and $R$ via a sandwich variance, and a Monte Carlo study.
The headline claims are: under correct specification, GLS $\le$ ANCOVA
$\le$ averaged in variance; the ANCOVA efficiency loss versus GLS is below
5% across a broad design regime; and under misspecification the simpler
estimators can outperform GLS.

## 2. Overall assessment and recommendation

**Recommendation: Reject in present form (major revision would be required
to reach even a resubmission, and the central empirical claims do not
survive verification).**

The manuscript is not ready for review at a top-tier statistics journal,
for three independent reasons, the first two of which are dispositive.

First, two of the three headline claims are contradicted by the paper's own
code. The claim that 'the efficiency loss from ANCOVA relative to GLS is
below 5% across a broad regime' is false: re-running the paper's own RE
table shows ANCOVA variance is typically two to ten times the GLS variance
(RE in the 0.1 to 0.5 range), with the 'below 5% loss' condition (RE $\ge$
0.95) holding in only about 5% of the tabulated grid, and essentially only
in the degenerate $J_1=1$ cases. The claim that the simpler estimators
'can outperform GLS' under misspecification is not demonstrated by either
misspecification figure; in both, misspecified GLS remains roughly an order
of magnitude more efficient than the averaged estimator at every parameter
value plotted.

Second, the manuscript is a working draft, not a finished paper. The
Simulation study section describes a 5000-replicate factorial Monte Carlo
producing empirical variance, power, and coverage, and the Reproducibility
section reports a specific master seed and per-replicate seeding scheme. No
such simulation exists in the source. The simulation chunk defines three
functions and runs a single trial under `set.seed(42)`. The SUR section
(4 lines plus a TODO), the 'When does the gap exceed 5%' section (a TODO
only), the 'When to use which estimator' section (a TODO), and the entire
Discussion (a TODO comment block) are unwritten.

Third, even fully executed, the novelty is thin relative to Frison and
Pocock (1992), Senn (2006), and the run-in-specific mixed-model work of
Wang et al. (2019) and Hu, Mackey and Thomas (2021); see Section 3.

I detail the verifiable problems below because they bear on the validity of
the project as a whole, including Paper 1, on which the GLS baseline rests.

## 3. Significance and novelty (skeptical)

The POST vs CHANGE vs ANCOVA efficiency ordering with one or more
pre-treatment measurements is the central result of Frison and Pocock
(1992, Statistics in Medicine), which establishes ANCOVA as the method of
choice with smaller variance than both alternatives and gives sample-size
formulas. Senn (2006) and Van Breukelen (2006) are the standard modern
treatments. The GLS / random-slopes power machinery for run-in trials is
covered by Wang et al. (2019) and by Hu, Mackey and Thomas (2021) (the
latter is in the bibliography). The only genuinely new element is the
specific closed-form RE of the three estimators expressed in the run-in
design parameters $(J_0, J_1, J_2)$ and $\sigma_b/\sigma$. That is an
incremental, special-case contribution. It would be of modest interest to
an applied trials-methodology readership (e.g., a short methods paper) if
the numbers were correct and the comparison complete, but it does not meet
the originality bar of Biometrika, and as written it falls short of
Statistics in Medicine as well. Critically, the one differentiating claim
that would have lifted significance, that ANCOVA is nearly as good as GLS
('below 5% loss'), is the claim that fails verification.

Searches performed (reported per instructions):
- WebSearch: 'Frison Pocock 1992 repeated measures mean summary statistics
  ANCOVA change score efficiency clinical trials'. Confirms the POST/CHANGE/
  ANCOVA ordering and the multi-baseline design result predate this paper.
- WebSearch: 'GLS REML efficiency loss covariance misspecification
  longitudinal robust sandwich variance Liang Zeger'. Confirms the standard
  framing (Liang and Zeger 1986; GEE robust sandwich) for the
  misspecification analysis the paper attempts.
- WebSearch: 'run-in baseline observations clinical trial GLS mixed model
  relative efficiency ANCOVA random slopes rate of change Alzheimer power'.
  Surfaces Wang et al. (2019) two-period LMM for run-in and Hu/Mackey/Thomas
  (2021), confirming the run-in GLS power story is established.

## 4. Major comments (correctness first; problem, why, remedy)

**M1. The abstract's central quantitative claim is false by the paper's own
code (Abstract lines 35-39; Results table chunk `re_table`, lines 336-401).**
The abstract states the ANCOVA-vs-GLS efficiency loss is 'below 5% across a
broad regime of $(J_0,J_1,J_2)$ and signal-to-noise ratios.' I reproduced
the paper's RE table (`re_ancova_vs_gls`) over the tabulated grid
($r\in\{0,1,2,4\}$, $k\in\{1,2,3\}$, $f\in\{0,1\}$, $\sigma_b/\sigma \in
\{0.5,1,2,5,20\}$). VERIFIED result: RE(GLS/ANCOVA) ranges from 0.000 to
1.000; it equals 1.000 only in the degenerate $(1,1,0)$ cell; for the
representative $(2,2,0)$ design RE is 0.32 to 0.36, i.e. ANCOVA variance is
roughly three times the GLS variance; the 'loss below 5%' condition (RE
$\ge 0.95$) holds in 5.2% of the grid (and only where $J_1=1$). Why this
happens (INSPECTED and VERIFIED by Monte Carlo, ratio empirical/formula =
1.007 at $(2,2,0)$): `var_gamma_ancova` collapses the $k$
post-randomisation visits into a single scalar mean before regression,
discarding exactly the within-period slope information that GLS exploits.
The ANCOVA-on-means estimator is therefore structurally far less efficient
than GLS whenever $J_1 > 1$, which is the regime of interest. The formula
is correctly implemented; the verbal claim about it is wrong. Remedy: delete
the 'below 5%' claim, or replace the ANCOVA-on-means estimator with an
ANCOVA that retains all post visits (e.g., per-visit response with run-in
covariate, or the cLDA estimator), for which a near-GLS efficiency claim
might actually hold. As written, the abstract, the Results section, and the
Conclusions ('ANCOVA is the recommended default') are unsupported.

**M2. The robustness claim is unsupported by the figures that purport to
support it (Abstract lines 39-44; Section 'Robustness to misspecification',
chunks `misspec_G` lines 418-491 and `misspec_R` lines 501-581).** The
abstract claims that under misspecification 'the simpler estimators can
outperform GLS' and 'ANCOVA and averaged remain comparatively robust.'
VERIFIED by re-running both chunks at the plotted settings ($r=2,k=2,f=0$,
$a=b=t=1$): in `misspec_G`, sweeping the assumed $\sigma_b^2$ over
$[0.01,5]$ with true $\sigma_b^2=1$, the misspecified-GLS sandwich variance
peaks at 5.30, while the averaged-estimator variance is 20.0; misspecified
GLS never exceeds averaged. In `misspec_R` (compound symmetry assumed,
AR(1) true), GLS sandwich variance ranges 0.74 to 2.16 across
$\rho\in[0,0.9]$, while averaged ranges 18.9 to 20.5; GLS is roughly
ten-fold more efficient at every $\rho$. The figures thus demonstrate the
opposite of the claim: GLS retains a large efficiency margin even when
badly misspecified, because the simpler estimators are so inefficient to
begin with. Why the claim is structurally hard to support here: for a
simpler estimator to overtake GLS, the GLS efficiency advantage (which M1
shows is large) must be erased by misspecification bias in the GLS variance;
the chosen misspecifications are far too mild to do that. Remedy: either
retract the robustness claim, or construct a misspecification regime where
it actually occurs (and note that the relevant comparison is misspecified
GLS variance against a correctly analysed simpler estimator, with the t-test
small-sample behaviour included), and report the crossover explicitly.

**M3. The Monte Carlo study does not exist (Section 'Simulation study',
lines 584-702; Reproducibility, lines 817-820).** The text states 'For each
configuration, 5000 trials are simulated. The outcomes are empirical
variance, empirical power, and coverage of nominal 95% confidence
intervals,' over a four-by-three-by-three-by-two factorial design. INSPECTED:
the `simulation` chunk defines `simulate_trial`, `analyze_gls`, and
`analyze_avg`, then runs exactly one trial under `set.seed(42)` and prints a
single point estimate and SE. There is no replicate loop, no factorial
sweep, no variance/power/coverage computation, and no output table. The
misspecified-GLS and ANCOVA analysers are never invoked in a simulation.
Consequently none of the described simulation results appear in the
manuscript. Separately, the Reproducibility section states the simulation
'sets `set.seed(20260418)` at its entry point' and derives 'per-replicate
seeds ... by `+ rep_idx`'; VERIFIED that neither `20260418` nor `rep_idx`
appears anywhere in the source. This is a fabricated reproducibility
statement describing code that was never written. This must be corrected
before any submission: either implement the simulation as described and
report results with Monte Carlo standard errors, or remove the simulation
section and the corresponding methods/abstract sentences entirely. As it
stands the section misrepresents work that does not exist.

**M4. No Monte Carlo standard errors, and no validation of the analytic
formulas in the manuscript itself.** Even the single validation trial
reports only a point estimate. There is no quantification of MC error
anywhere (necessarily, given M3). The analytic RE claims, which are the
core of the paper, are presented without the simulation cross-check the
paper advertises. (For the record, I independently verified the GLS,
ANCOVA, and averaged variance formulas by direct simulation: GLS empirical
0.0221 vs formula/m 0.0216 at $m=100$; ANCOVA empirical 0.0336 vs formula/m
0.0333 at $m=200$; the Frost closed form $(\sigma^2+2\sigma_b^2t^2)/t^2$
matches `var_gamma_matrix(0,2,0)` exactly. The formulas are sound; the
manuscript simply does not demonstrate this.) Remedy: report MC SEs for
every simulated quantity and include a formula-vs-simulation validation
table.

**M5. The MIRIAD numerical example produces non-credible, uncommented
sample sizes (chunk `miriad_comparison`, lines 709-770).** With
$\sigma^2=0.0025$, $\sigma_b^2=0.9745$ (so $\sigma_b/\sigma \approx 19.7$),
$\Delta=0.25$, 90% power, VERIFIED that the table yields GLS $N=1$ per group
for many designs (e.g. $(4,3,1)$, $(2,3,0)$) alongside averaged $N$ up to
8193. An $N=1$-per-group 'sample size' is not a usable design recommendation
and signals the example was not sanity-checked; it arises because at this
extreme SNR the per-pair GLS variance is minuscule, but the table presents
it without comment next to four-figure averaged sample sizes. This same
table also flatly contradicts M1's 'below 5%' claim: ANCOVA needs $N$ up to
16 against GLS $N=1$. Remedy: justify the MIRIAD variance components and
$\Delta$ on the original (non-standardised) scale, impose a floor that makes
$N$ meaningful, and remove or reconcile the example with the efficiency
claims.

**M6. The ordering claim is correct but under-stated; the proof is asserted,
not given (Section 'When does the gap exceed 5%', and BLUP section lines
282-319).** VERIFIED numerically across 413 grid points: GLS $\le$ ANCOVA
$\le$ averaged in variance holds with zero violations. This is the one
headline claim that survives. However, the paper only gestures at the
argument (Gauss-Markov optimality of GLS; the averaged and ANCOVA
estimators as linear unbiased estimators in the same model). A top-tier
paper should state and prove this as a proposition: both simpler estimators
are linear unbiased estimators of $\gamma$ under the random-slopes model, so
each has variance at least that of the GLS/BLUE; and ANCOVA $\le$ averaged
because ANCOVA additionally optimises over the run-in coefficient $\phi$
while averaged fixes it at the change-score contrast. The 'When does the gap
exceed 5%' section that should contain this is an empty TODO (line 405).
Remedy: supply the proposition and proof, and characterise the boundary
analytically rather than deferring it.

**M7. Unfinished sections presented as part of a complete manuscript.** The
SUR subsection (lines 259-279) introduces a seemingly-unrelated-regression
estimator, states a reparameterisation, then stops at a TODO to 'derive the
SUR variance and compare'; SUR is nonetheless referenced in the Discussion
TODO as 'a middle ground.' The Discussion (lines 779-790) is entirely a
comment block. 'When to use which estimator' (lines 772-776) is a TODO.
These are not minor: roughly a third of the paper's intended content is
absent. Remedy: complete or remove. If SUR is not analysed, remove the
subsection rather than leaving a stub in a submitted manuscript.

**M8. cLDA is named as the more relevant competitor but excluded (lines
136-149).** The paper correctly notes that the constrained longitudinal data
analysis estimator (Liang and Zeger 2000; Lu 2010) sits between averaged and
GLS, exploits the full covariance without specifying $G$, and is the
natural robust alternative, then declares a full comparison 'beyond the
scope.' Given that M1 and M2 show ANCOVA-on-means and averaged are both
poor, and that cLDA is precisely the estimator that would plausibly deliver
the near-GLS-with-robustness story the abstract wants to tell, excluding it
removes the paper's strongest potential result. Remedy: include cLDA; it is
arguably the point of the comparison.

## 5. Minor comments

- m1. Notation drift: the design code uses $(r,k,f)$ throughout while the
  prose uses $(J_0,J_1,J_2)$. State the mapping ($J_0=r$, $J_1=k$, $J_2=f$)
  once, explicitly, and use one convention in the equations.
- m2. Equation (avg, lines 112-117) defines $d_{ik}$ with per-subject limits
  $J_{2,i}, J_{0,i}$ implying subject-varying visit counts, but the entire
  development and code assume common, fixed $J_0,J_1,J_2$. Remove the
  subject indices or address unbalanced designs explicitly.
- m3. The ANCOVA variance formula uses $(1-\rho^2)$ with $\rho$ the marginal
  correlation of run-in and post means. State that this is the population
  ANCOVA variance with $\phi$ known; the feasible ANCOVA estimates $\phi$
  and carries the usual finite-sample inflation, which should be acknowledged
  alongside the analogous feasible-GLS caveat in the Introduction.
- m4. Reproducibility section claims `nlme` is 'the GLS fitter for the
  simulation comparator.' VERIFIED `nlme` is never called; `analyze_gls`
  implements GLS by hand. Correct the statement.
- m5. Reproducibility section cites `tinytest` as the testing harness; the
  project's CLAUDE.md and Makefile reference `testthat`. Reconcile.
- m6. Data availability and Reproducibility refer to paths
  `analysis/paper{1,2,3}-*/` and `analysis/paper4-gls-vs-ttest/`; the actual
  layout is `analysis/report/0{1,2,3,4}-*/`. Fix the paths.
- m7. Abstract calls GLS 'the minimum-variance unbiased estimator when
  $\Sigma$ is known'; it is the minimum-variance estimator among linear
  unbiased estimators (BLUE), not MVUE in general. Tighten.
- m8. Figure captions for `misspec_G`/`misspec_R` assert the averaged
  estimator is 'unaffected' by misspecification; it is unaffected in
  estimator definition but its variance does change with the true $R$
  (visible in `v_avg_ar1`). Reword to 'its variance does not depend on the
  assumed model.'
- m9. The single-trial validation output (`set.seed(42)`) is left in the
  rendered document; remove debugging output from the final manuscript.

## 6. Missing or questionable references

| Key / item | Issue | Action |
|---|---|---|
| `frison1992` | In bib, the foundational POST/CHANGE/ANCOVA efficiency result with multiple baselines, NOT cited in this paper despite being the closest prior work | Cite prominently in Introduction and positioning |
| `frison1997` | In bib, linearly divergent (slope) treatment effects via summary statistics, directly relevant to the slope model here, not cited | Cite |
| `hu2021` | Defined THREE times in `references.bib` with conflicting metadata (Biometrical Journal 63(4) vs Int J Biostatistics; author 'Nan'/'Na' Hu; 'Howard'/'Hugh' Mackey). BibTeX will silently take one and warn | De-duplicate into two distinct keys (the BJ monotone-missing paper and the IJB power-formula paper) |
| `liang1986` | The canonical robust-sandwich / misspecification reference underlying Section 4, not cited where the sandwich variance is introduced | Cite at the sandwich construction |
| `lu2009` (Lu, Mehrotra, Liu, sample size for cLDA) | In bib, relevant to the cLDA discussion, not cited | Cite if cLDA is added (M8) |
| `crager1987` | Cited (line 128) but the claim it supports (formal efficiency comparison) is more precisely Frison-Pocock; verify attribution | Check |
| `senn1997`, `senn2002`, `jones2014`, many AD-imaging entries | Present in the shared 60-entry bib but unused by Paper 4 | Acceptable for a shared bib, but if the bib is paper-specific at submission, prune to cited items |

I did not find any cited reference that fails to exist; the problems are
miscitation (M-level), omission of the closest prior work (Frison-Pocock),
and duplicate keys.

## 7. Suggestions

- Reframe as what the verified evidence supports: under correct
  specification GLS strictly dominates, and ANCOVA-on-means and averaged
  carry substantial efficiency penalties that grow with $J_1$; quantify the
  penalty rather than claiming it is negligible.
- If the intended message is 'simpler estimators are robust,' build it on
  cLDA (M8) and on a finite-sample comparison that includes t-test
  small-sample inflation, not on the inefficient averaged estimator.
- Implement the advertised simulation with MC SEs and a formula-validation
  table; the validation I ran shows the formulas are correct, so this is
  achievable and would materially strengthen the paper.
- State the ordering result as a proved proposition (M6).
- Resolve the duplicate `hu2021` keys and add Frison-Pocock before any
  resubmission.

## 8. Scope of review

I read the full .Rmd including all code chunks and the shared bibliography,
read the relevant package source it depends on, and read the Mathematica
derivations (`woodbury_variance.m`, `blup_derivation.m`). I re-executed the
paper's variance functions, RE tables, and both misspecification chunks
under R 4.5.3 via `devtools::load_all`, and independently validated the
GLS, ANCOVA, and averaged variance formulas by direct Monte Carlo
(claims tagged VERIFIED above rest on these runs; INSPECTED claims rest on
reading the source). I did not attempt a full PDF render, did not assess
the CSL/formatting, and did not independently re-derive the Woodbury closed
forms by hand (these belong to Paper 1; I confirmed numerical agreement
with the package and with the Frost closed form). The novelty assessment
rests on three WebSearch queries (Section 3) plus the cited literature; I
did not exhaustively survey the run-in design literature.

---
*Rendered on 2026-06-13 at 18:18 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/analysis/report/04-gls-vs-ttest/referee-report-2026-06-13.md*
