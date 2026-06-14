# Referee Report: 'Power Analysis for Run-In Clinical Trials under AR(1) and General Correlation Structures'
*2026-06-13*

Manuscript: `analysis/report/02-correlation/report.Rmd` (Paper 2 of a
four-paper compendium). Author: Ronald G. Thomas (UCSD).

Confidence labels used throughout: VERIFIED (I executed code or a symbolic
check and observed the result), INSPECTED (I read the source and confirmed
by reasoning), INFERRED (consistent with surrounding material but not
directly executed), SUSPECTED (a probable defect I could not fully
isolate).

---

## 1. Summary

The manuscript proposes to extend the closed-form run-in power results of
Paper 1 (compound symmetry, CS) to AR(1) and more general correlation
structures. It states two contributions: (i) a Woodbury-identity
derivation of the variance of the treatment-effect estimator under an
AR(1) marginal covariance, accommodating arbitrary numbers of run-in
($J_0$) and common-close ($J_2$) observations; and (ii) a moment-based
derivation via the cross-covariance between run-in and treatment-period
means. The abstract asserts that the two routes 'yield identical results,
providing internal validation', develops a 'taxonomy of five correlation
cases', and reports that AR(1)-corrected sample sizes exceed the
compound-symmetric ones by 5 to 25 percent, with the largest deviation at
small $J_0$ and small $\rho$.

I checked every quantitative claim against the manuscript's own code and
against independent re-derivations in R. The Woodbury (GLS) implementation
is correct. The remaining headline claims are not supported by the
manuscript as written: the moment-based variance formula is wrong, the
two-route equivalence is never actually demonstrated and the only numerical
comparison in the paper shows the two routes disagreeing, and the central
5 to 25 percent sample-size claim is contradicted in both direction and
magnitude by running the paper's own functions. The manuscript also
contains four visible `TODO` placeholders in core sections, including the
derivation that the abstract describes as completed.

## 2. Overall assessment and recommendation

**Recommendation: Reject in present form (major revision required before
any reconsideration).**

The paper is a working draft, not a finished manuscript. Three of the four
results sections contain unfinished derivations marked with `TODO`
comments (Sections 2.2, 6, 8), the AR(1) change-score correlation
subsection literally reads '(to be derived)', and the Discussion is a
single `TODO` line. Independently of the incompleteness, the two central
scientific claims (two-route internal validation; the 5 to 25 percent
correction) are demonstrably false as stated. The novelty premise is also
overstated: the paper claims the standard CRAN tool provides no analytic
AR(1) formula, but `longpower::power.mmrm.ar1` is exactly such a formula.

I would be willing to look at a substantially rewritten version in which:
the moment-based formula is corrected and the two estimands are
reconciled; the equivalence claim is either proven for a common estimand
or withdrawn; the sample-size claim is recomputed and the direction of the
AR(1) effect is correctly characterised; the placeholders are completed;
and the positioning against `longpower`, Frison-Pocock, Winkens et al. and
Munoz et al. is corrected. That is a major-revision-to-resubmission
amount of work.

This manuscript is not at the top-tier bar (Biometrika / Statistics in
Medicine methods). On a corrected and completed basis the contribution
would be a competent applied-methods note; tier fit would be a specialist
trials journal rather than a flagship.

## 3. Significance and novelty

The applied question (does the CS-based run-in efficiency gain survive when
the truth is AR(1)?) is reasonable and of practical interest for AD
neuroimaging trials. However the novelty as framed is thin and partly
misstated.

- The claim in the Introduction (lines 114-127) that `longpower` provides
  'No analytic AR(1)-residual formula' is incorrect. The package exports
  `power.mmrm.ar1`, an analytic (non-simulation) sample-size routine for a
  mixed model with AR(1) correlation, based on Lu, Luo and Chen (2008).
  [VERIFIED via package documentation and web search.] The cited
  `iddi2022` R Journal paper documents this function. The stated gap
  therefore does not exist in the form claimed. The genuine increment, if
  any, is the combination of AR(1) residuals with the run-in / common-close
  augmented design and the random-slope structure, which `power.mmrm.ar1`
  does not directly address. The paper must reframe around that narrower
  and defensible claim.

- The moment-based route (Section 4) is, structurally, the Frison-Pocock
  (1992) summary-statistics approach specialised to AR(1). Frison-Pocock is
  in the bibliography (`frison1992`) but is not cited in Section 4, and the
  paper does not position its phase-mean variance derivation relative to
  that established work. As written, Section 4 reads as a rediscovery.

- Winkens et al. (2007, `winkens2007`, cited) already established the key
  qualitative result this paper rediscovers numerically: under AR(1) the
  marginal value of additional repeated measures is small, whereas under CS
  it is large. The paper should engage this rather than present the
  robustness figure as new insight.

Verdict: incremental and, in its current framing, partly non-novel.
Genuine residual contribution exists (run-in plus AR(1) plus random slope
in closed form) but is neither isolated nor correctly positioned.

## 4. Major comments (correctness first)

**M1. The two-route 'internal validation' is never demonstrated and the
paper's own numbers contradict it. [VERIFIED]**
Problem. The abstract (lines 22-25) and Methods state 'Both routes yield
identical results, providing internal validation.' This is the headline
argument of the paper. It is not shown anywhere, and the single numerical
comparison that exists (chunk `moment_based`, lines 295-326) prints two
numbers that disagree: moment-based 1.4375 versus GLS 1.6662. [VERIFIED by
executing the chunk's code.] The chunk's own caption excuses the mismatch
('Moment-based ignores random slope'), which concedes that the two routes
do not compute the same quantity.
Why it matters. The two routes target different estimands: `var_gamma_ar1`
returns $\Var(\hat\gamma)$ for the slope-difference fixed effect $\gamma$
(scaled by $t$) under the full random-slopes GLS model, whereas
`ar1_var_change` returns the variance of a phase-mean difference
$\bar X_K - \bar X_0$ with no random slope. These are not the same number
and cannot validate each other. Section 5 (lines 329-353) in fact states
plainly that 'The GLS approach ... and the moment-based approach ...
compute different quantities', directly contradicting the abstract.
Remedy. Either (a) define a common estimand (for example set
$\sigma_b^2 = 0$, where GLS should reduce to a generalised-least-squares
weighting of the same change scores, and prove algebraic equality of the
two closed forms), or (b) withdraw the equivalence claim entirely and
reframe the moment-based route as an alternative estimator with its own
variance. Note that the natural check in (a) cannot currently be run
because `var_gamma_ar1` throws on $\sigma_b^2 = 0$ (see M3).

**M2. The moment-based variance-of-phase-mean formula is incorrect.
[VERIFIED]**
Problem. The closed form for $S_n$ (lines 291-293) and its use in
`ar1_var_mean` (lines 303-308) are inconsistent. The stated definition is
$S_n = \sum_{i=1}^{n-1}(n-i)c^{i-1}$, and the closed form
$S_n = [n(1-c)-(1-c^n)]/(1-c)^2$ correctly evaluates that sum. [VERIFIED.]
But the variance of the mean of $n$ equally spaced AR(1) observations is
$\Var(\bar X_n) = \frac{\sigma^2}{n^2}\big(n + 2\sum_{l=1}^{n-1}(n-l)c^{l}\big)$,
which requires $c^{l}$, not $c^{l-1}$. The implemented expression
$\frac{\sigma^2}{n}\big(1 + \frac{2\rho}{nc}S_n\big)$ reduces to
$\frac{\sigma^2}{n^2}\big(n + 2\sum (n-i)c^{i-1}\big)$, i.e. it is off by
exactly the factor $c$ inside the lag sum.
Why it matters. The error is large and not cosmetic. For $n=2,\ \rho=0.3$
the formula returns $1.000\sigma^2$ versus the correct $0.650\sigma^2$; for
$n=4,\ \rho=0.5$ it returns $0.781\sigma^2$ versus $0.516\sigma^2$.
[VERIFIED against brute-force $\mathbf 1'\Sigma\mathbf 1/n^2$.] Because the
cross-covariance term `ar1_cov_phase_means` is correct (it matches
brute-force exactly, VERIFIED), the change-score variance `ar1_var_change`
inherits the variance-term error and is wrong by a similarly large margin.
Remedy. Replace $c^{i-1}$ by $c^{i}$ in the lag sum (equivalently drop the
spurious $1/c$ factor): use
$\Var(\bar X_n) = \frac{\sigma^2}{n^2}\big(n + 2\,T_n\big)$ with
$T_n = \sum_{l=1}^{n-1}(n-l)c^{l} = \frac{c\,[\,n(1-c)-(1-c^n)\,]}{(1-c)^2}$.
Re-derive Section 4 and re-run all dependent output.

**M3. `var_gamma_ar1` fails at $\sigma_b^2 = 0$, blocking the only valid
equivalence check. [VERIFIED]**
Problem. `var_gamma_ar1` forms $H = (G^{-1} + Z'R^{-1}Z)^{-1}$ via
`solve(solve(G))` with $G = \sigma_b^2$. At $\sigma_b^2 = 0$ this throws
'system is exactly singular'. [VERIFIED.] The natural limit
($\sigma_b^2 \to 0$, where GLS reduces to ordinary change-score GLS and
should match a correctly derived moment route) is therefore not
computable.
Why it matters. The boundary case is exactly the one Section 5 invokes to
argue the two routes 'coincide'. The claim cannot even be tested with the
current code.
Remedy. Use the algebraically equivalent form
$H = G - G Z'(R + ZGZ')^{-1} Z G$, or guard the $\sigma_b^2 = 0$ case
analytically, so the limit is well defined.

**M4. The central 5 to 25 percent sample-size claim is contradicted by the
manuscript's own functions, in both direction and magnitude. [VERIFIED]**
Problem. The abstract (lines 32-36) and Conclusions state AR(1)-corrected
sample sizes exceed CS ones by 5 to 25 percent, largest at small $J_0$ and
small $\rho$. Executing the paper's `comparison` chunk (lines 405-437) with
the stated MIRIAD parameters ($\sigma^2=0.0025$, $\sigma_b^2=0.9745$,
$t=1$, $\Delta=0.25$, $\rho=0.5$, 90 percent power) yields, for
$J_0 = 0,1,2,3,4$: $N_\text{CS} = 329, 3, 2, 1, 1$ and
$N_\text{AR(1)} = 328, 2, 1, 1, 1$. [VERIFIED.] The AR(1) figures are
lower than or equal to the CS figures, not 5 to 25 percent higher. More
fundamentally, the variance ratio $V_{AR(1)}/V_{CS}$ is below 1 for
essentially all $\rho > 0$ (e.g. $J_0=2$: $0.77$ at $\rho=0.5$, $0.18$ at
$\rho=0.9$; VERIFIED), because the implemented AR(1) change-score
covariance coincides with the CS covariance at $\rho = 0$ and the variance
decreases as $\rho$ grows. The qualifier 'largest at small $\rho$' is
therefore also backwards: deviation grows with $\rho$.
Why it matters. This is the paper's only quantitative conclusion and it is
wrong as stated. The direction of the AR(1) correction is the opposite of
what the abstract claims under the implemented model.
Remedy. Recompute and re-state. Decide whether the intended comparison
holds $\rho$ fixed across structures or matches a common lag-1 correlation
(the two give very different stories), and characterise the sign of the
effect correctly. Separately, the MIRIAD sample sizes of 1 to 3 are
implausible and indicate that $\Delta = 0.25$ is far too large for this
variance scale; reconcile units (this issue is shared with Paper 1 and
should be fixed compendium-wide).

**M5. AR(1) marginal covariance: formula is correct but the limit claim is
conceptually wrong. [VERIFIED / INSPECTED]**
The change-score covariance
$R_{jl} = \sigma^2(\rho^{|j-l|} - \rho^{|j|} - \rho^{|l|} + 1)$ (lines
186-189, implemented in `build_R_matrix_ar1`) matches the brute-force
contrast $L\,\Sigma_\text{raw}\,L'$ exactly. [VERIFIED.] However the
Introduction (lines 124-127) claims consistency with `longpower` 'when the
correlation is replaced by its compound-symmetry limit'. CS is not a limit
of AR(1): the implemented model collapses to CS at $\rho \to 0$
(VERIFIED: `build_R_matrix_ar1(..., 1e-9)` equals
`build_R_matrix(...)`), which is the no-serial-correlation case, not a
compound-symmetric limit. The sentence conflates two distinct structures
and should be corrected.

**M6. Section 2.2 change-score formulas silently switch the change-score
definition. [VERIFIED]**
Problem. Section 2.2 declares the change-from-baseline parameterisation
$C_j = Y_j - Y_0$ (line 159), then gives
$\Corr(C_1,C_3) = bt^2/(2a+bt^2)$ and
$\Corr(C_1,C_2) = (-a+bt^2)/(2a+bt^2)$ (lines 166-170). Under
change-from-baseline these are wrong, but under the consecutive-increment
parameterisation $C_j = Y_j - Y_{j-1}$ they are exactly correct.
[VERIFIED both ways.] The AR(1) machinery elsewhere uses
change-from-baseline. The section therefore mixes two incompatible
definitions without comment, and never defines the time points to which
$C_1, C_2, C_3$ correspond.
Remedy. Fix the parameterisation consistently, state the model and the
times explicitly, and reconcile with the change-from-baseline convention
used in `build_R_matrix_ar1`.

**M7. Core sections are unfinished. [INSPECTED]**
The AR(1) change-score correlation derivation reads '(to be derived)' with
a `TODO` (lines 172-175); the stratified-variance section ends in a `TODO`
(lines 397-398); the Discussion is a single `TODO` (lines 475-476). The
abstract nonetheless describes the derivation as completed and the taxonomy
as 'developed'. A submission cannot claim completed results that the body
flags as pending.

**M8. The 'five-case taxonomy' is not load-bearing. [INSPECTED]**
The taxonomy (lines 132-155) lists CS, AR(1), CS-with-variable-run-in,
AR(1)-with-variable-run-in-and-close, and unstructured. Cases 3 and 4 are
just Cases 1 and 2 with subject-varying $J_0, J_2$ (a stratification, not a
correlation structure), and Case 5 (unstructured) is never derived or used.
The abstract promises the taxonomy 'characterises regimes in which run-in
efficiency gains are robust', but no case boundaries or regime thresholds
are derived anywhere. Either derive the boundaries (e.g. the $\rho$ and
$\sigma_b/\sigma$ regions where the CS and AR(1) variances agree to within
a tolerance) or drop the taxonomy framing.

## 5. Minor comments

- The GLS implementation `var_gamma_ar1` is correct: it matches an
  independent block-diagonal GLS computation across $J_0 \in \{0,2,4\}$ and
  $\rho \in \{0.3, 0.6\}$ to machine precision. [VERIFIED.] This should be
  stated as the paper's one solid result.
- Data-availability and Reproducibility sections reference paths that do
  not match the repository layout: they cite `analysis/paper2-correlation/`
  and `analysis/paper{1,3,4}-*/`, but the actual directories are
  `analysis/report/02-correlation/` etc. Fix before submission.
- `references.bib` is a 60-plus-entry shared collection imported for Zotero,
  not a curated per-paper list; only seven keys are actually cited
  (`frost2008`, `wang2019`, `winkens2007`, `diggle2002`, `verbeke2000`,
  `chambless1993`, `iddi2022`). [VERIFIED by grep.] The bib contains
  duplicate `hu2021` keys (three competing definitions, lines 255, 505,
  701) which will cause a BibTeX clash. Trim to a per-paper bibliography.
- The `moment_based` chunk passes `ar1_var_change(2, 2, 1, a, 0.5)` with
  `a = 1.0` as $\sigma^2$ and `0.5` as $\rho$; this is correct positionally
  but fragile (no named arguments). Use named arguments throughout.
- Reproducibility section claims `set.seed(20260418)` is set 'at the
  driver's entry point' for Monte Carlo benchmarking, but no Monte Carlo
  benchmarking appears in the manuscript. Remove or implement.
- `geometry: left=3cm, right=5cm` produces a 5 cm right margin, leaving an
  unusually narrow text block. Likely an unintended copy from a template.
- Notation: the moment-based section introduces $c = \rho$ (line 271) while
  Paper 1 and the Mathematica sources use $c = \sigma^2 - \sigma_b^2 t^2$
  as a reparameterisation. Reusing $c$ for the autocorrelation invites
  confusion across the compendium. Choose a distinct symbol.

## 6. Missing or questionable references

| Reference | Status | Issue / required action |
|---|---|---|
| Lu, Luo, Chen (2008), `lu2008` | In bib, not engaged on AR(1) | This is the basis of `longpower::power.mmrm.ar1`, the analytic AR(1) routine the paper claims does not exist. Must be discussed and contrasted. |
| Frison & Pocock (1992), `frison1992` | In bib, not cited | Canonical summary-statistics / phase-mean variance method and the 'more than one pre-treatment measurement' result. Section 4 must cite and position against it. |
| Munoz, Carey, Schouten, Segal, Rosner (1992), Biometrics 48:733-742 | Missing | Foundational parametric serial-correlation family for longitudinal data; standard citation when motivating AR(1) over CS. Add. |
| Winkens et al. (2007), `winkens2007` | Cited once | Already established that AR(1) reduces the value of intermediate measures; the robustness figure rediscovers this. Engage substantively, not in passing. |
| Diggle, Heagerty, Liang, Zeger (2002), `diggle2002` | Cited generically | Should anchor the specific AR(1) marginal-covariance and serial-correlation discussion, not just a general pointer. |
| longpower `power.mmrm.ar1` (documented in `iddi2022`) | Mischaracterised | The Introduction's gap claim is false; correct it. |

Bib hygiene: duplicate `hu2021` keys; one `winkens2005` entry has DOI
`10.1002/sim.2370` while the web record for that title shows `10.1002/
sim.2385` (verify the correct DOI). The shared bib should be reduced to the
seven cited entries plus the additions above.

## 7. Suggestions (constructive, beyond mandatory fixes)

- Lead with the one correct and genuinely useful result: a verified
  closed-form AR(1) GLS variance for the run-in / common-close
  random-slope design. Build the paper around proving and characterising
  that, rather than the unsubstantiated two-route equivalence.
- For the comparison to be interpretable, fix the cross-structure matching
  convention. Comparing CS and AR(1) at a common nominal $\rho$ confounds
  marginal variance with correlation decay. Match on a common lag-1
  correlation or a common total information, and state which.
- Replace the five-case 'taxonomy' with an explicit derivation of the
  $(\rho, \sigma_b/\sigma, J_0)$ region in which the CS formula is within,
  say, 5 percent of the AR(1) formula. That would deliver the practical
  'map' the abstract promises.
- Add a small simulation that fits the assumed mixed model to AR(1)-
  generated data and recovers the analytic variance; this is the honest
  form of the 'internal validation' the paper currently claims.
- Resolve the units problem so the MIRIAD example yields realistic
  per-group sample sizes (tens to hundreds), and keep that resolution
  consistent across Papers 1 to 4.

## 8. Scope of review

I read the full manuscript including all eight sections, the YAML header,
and every R chunk. I read the shared bibliography
(`analysis/references.bib`, to which the paper's `references.bib` is a
symlink) and the available Mathematica derivation sources
(`analysis/scripts/woodbury_variance.m`, `raw_outcome_model.m`); the
`frost311.m` / `frost312.m` files named in the review brief do not exist in
the repository, and `docs/frost.pdf` is the Frost et al. (2008) source
paper rather than an internal derivation. I cross-read Paper 1
(`analysis/report/01-runin-power/report.Rmd`) for notation and parameter
consistency.

Verification performed in R by sourcing the package functions directly
(`R/design.R`, `R/covariance.R`, `R/variance.R`, `R/covariance-ar1.R`,
`R/power.R`): (i) brute-force confirmation of the AR(1) change-score
covariance; (ii) block-diagonal GLS confirmation of `var_gamma_ar1`;
(iii) brute-force confirmation that the moment-based variance-of-mean
formula is wrong and the cross-covariance formula is correct; (iv)
re-execution of the `comparison` chunk to test the 5 to 25 percent claim;
(v) confirmation of the $\sigma_b^2 = 0$ failure; (vi) confirmation of the
Section 2.2 parameterisation ambiguity. Literature checks used WebSearch
and WebFetch (Winkens et al. 2007; longpower `power.mmrm.ar1`; Munoz et al.
1992; Frison-Pocock 1992) and are reported in Sections 3 and 6.

I did not render the full PDF, re-run the Mathematica notebooks
symbolically, or audit Papers 3 and 4. The symbolic two-route equivalence
(M1) was assessed by estimand analysis and numerical contradiction rather
than a full closed-form derivation, because the relevant boundary case is
not computable in the current code (M3); a definitive symbolic proof of
non-equivalence was therefore not produced, though the numerical evidence
of disagreement is unambiguous.

---
*Rendered on 2026-06-13 at 18:17 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/analysis/report/02-correlation/referee-report-2026-06-13.md*
