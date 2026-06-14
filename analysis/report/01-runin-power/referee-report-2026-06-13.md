# Referee Report: 'Power Analysis for Longitudinal Clinical Trials with Run-In and Common Close Observations'
*2026-06-13*

## 1. Summary

The manuscript derives closed-form and matrix-computable expressions for the
variance of the treatment-effect (rate-difference) estimator in a two-arm
longitudinal trial with a random-slopes linear mixed model, analysed under a
change-score (equivalently cLDA) parameterisation. The derivation applies the
Woodbury matrix identity to the GLS information matrix. Three nested cases are
treated: (i) the conventional design with no run-in, recovering Frost, Kenward
and Fox (2008); (ii) an explicit two-run-in case; and (iii) a general
configuration with J0 run-in, J1 treatment, and J2 common-close observations,
evaluated numerically. Numerical illustrations use parameters attributed to
Alzheimer's disease neuroimaging trials.

The core algebra that I was able to verify is correct, and the recovery of the
Frost (2008) and Ard-Edland (2011) boundary cases is genuine and numerically
exact. However, the manuscript contains a substantive correctness error in
Section 6 (the J0=1 claim), a misappropriation of the Frost MIRIAD variance
parameters that renders the entire Section 5.5 application numerically
meaningless and internally contradictory, an unsupported headline efficiency
claim in the abstract, and an overstated novelty position relative to Frost
(2008), Wang et al. (2019), and the Hu-Mackey-Thomas (2021) pair. The
common-close generalisation is the only clearly novel contribution, and it is
presented without a closed form, without a comparison to delayed-start designs,
and without validation.

## 2. Overall assessment and recommendation

**Recommendation: major revision** (closer to reject-and-resubmit than to
minor). The verified algebra is sound and the framework is internally coherent,
which keeps this above the rejection threshold. But the application section is
quantitatively broken, one stated theoretical conclusion is false, and the
novelty as currently framed does not clear the bar for a top-tier methodological
journal (Biometrika / Statistics in Medicine methodology). In its present form
the genuinely new content (the common-close result) would fit a applied-design
or software-companion venue. With the corrections below and a defensible novelty
claim centred on the common-close phase, a Statistics in Medicine submission is
plausible; Biometrika would require substantially more (asymptotics, optimality,
or a result that is not a direct Woodbury bookkeeping extension of Frost).

## 3. Significance and novelty

What is actually new, stated precisely:

- The **explicit closed form for the two-run-in case** (Eq. var_r2k1) is, as far
  as I can determine, not in Frost (2008), which gives only the one-run-in case
  (its Eq. 16). This is a minor, incremental closed form. VERIFIED that it
  matches the matrix computation.
- The **general (J0, J1, J2) Woodbury evaluation including a common-close phase**
  is the principal novel element. I could locate no prior closed-form
  random-slopes GLS variance that incorporates a post-treatment common-close
  phase framed this way. This is plausibly novel but is delivered only as a
  numerical recipe, not a formula, and is not positioned against the
  delayed-start / randomized-withdrawal literature that readers will treat as the
  natural comparator.

What is not new, and is overclaimed:

- The J0 run-in result is the natural generalisation of Frost (2008). Frost
  already establishes the framework, the change-score covariance, the GLS
  variance route, the run-in design, and the MIRIAD application; Frost (2008)
  Section 2.1 even gives the general LMM sample-size machinery for 'any design
  matrix and postulated Sigma'. The present manuscript's J0 generalisation is
  Woodbury bookkeeping on top of that. The abstract's framing ('have remained
  scattered ... and partially specified') overstates the gap.
- The recovery of the Edland and Diggle formulae (Section 3.3) is a re-derivation
  of known boundary results, correctly done but not novel.
- The author's own prior work (Hu, Mackey, Thomas 2021, both papers) already
  provides closed-form power for random slope-and-intercept rate comparisons with
  missing data. The novelty delta against the author's own results must be stated
  crisply and is not.

Contemporariness: the engagement with post-2015 work is thin. Nash, Morgan,
Frost and Mulick (2021, Stata Journal, slopepower, plus its 2024 erratum) is the
direct software realisation of the competing Frost approach and is not cited.
Iddi and Donohue (2022, longpower) is cited only in passing; it should be the
explicit numerical benchmark.

## 4. Major comments (correctness first)

**M1. Section 6, lines 934-944: the claim that Var(gamma-hat)|J0=1 can exceed
Var(gamma-hat)|J0=0 is FALSE for this estimator. [VERIFIED]**
The manuscript states that adding a single run-in observation introduces delta
and 'can offset the information gained ... so that Var|J0=1 may exceed Var|J0=0
for some parameter configurations', and concludes designers should plan zero or
at least two run-in observations, never one. I evaluated var_gamma_matrix(0,k,0,
.) versus var_gamma_matrix(1,k,0,.) exhaustively over k in 1..6, sigma_b^2 in
[1e-4, 1e3] (40 log-spaced points), and t in {0.25, 0.5, 1, 2, 4}. In no
configuration does J0=1 exceed J0=0. At sigma_b^2 = 0.5, t=1, k=2 the two are
exactly equal (both 2.000); elsewhere J0=1 is strictly smaller. The correct
statement is that adding a single run-in observation never increases the
variance but at certain parameter values yields no improvement (the run-in
information exactly cancels the cost of the delta degree of freedom). The
practical recommendation as written ('not one') is therefore unsupported.
*Why it matters:* this is presented as a design conclusion and attributed to
Frison-Pocock (1997); it would mislead trial designers. *Remedy:* either prove
the inequality the author believes holds (I believe it does not for the
change-score GLS estimator) or replace the paragraph with the correct
'never-worse, sometimes-no-gain' characterisation, with the boundary condition
under which the gain is zero. The Frison-Pocock citation should be dropped or
re-scoped, as their result concerns a different estimator and parameterisation.

**M2. Section 5.5, lines 1093-1156: the Alzheimer's application misuses the Frost
MIRIAD variance components, producing absurd and self-contradictory sample
sizes. [VERIFIED against frost.pdf Table II]**
The manuscript sets sigma^2 = 0.0025 (residual) and sigma_b^2 = 0.9745. Frost
(2008) Table II reports three components for MIRIAD: between-subject slope
variance sigma_b^2 = 0.9745; within-subject between-visit effects sigma_u^2 =
0.3338; and additional measurement error sigma^2 = 0.0025. The manuscript has
taken sigma_b^2 correctly but has used Frost's tiny 'additional measurement
error' (0.0025) as the change-score residual and has silently discarded
sigma_u^2 = 0.3338, which is the dominant within-subject term and is precisely
what plays the role of the residual in a change-score analysis. The consequence:
the AD example yields n = 329 for the standard design but n = 2 for (J0=2, J1=2,
J2=0) and n = 1 for the full design, i.e. 99 to 100 percent 'reductions' (lines
1148-1155). These are not plausible trial sizes and they directly contradict the
abstract's 30-45 percent claim and Frost's own Table III (which reports 175 per
arm for a no-run-in 12-month design, 122 for 24 months). With a defensible
residual of sigma_u^2 + sigma^2 = 0.3363, I obtain n = 385 (standard) and n =
145 (J0=2), a 62 percent reduction, which at least lies in a sane range.
*Why it matters:* the only empirical demonstration in the paper is wrong by two
orders of magnitude and undermines every quantitative claim. *Remedy:* reconcile
the change-score residual with Frost's variance decomposition. State explicitly
how sigma_u^2 maps into the change-score covariance R. If the differenced-error
covariance R = sigma^2(I + 11') is meant to absorb only measurement error, then
the model is omitting the visit-level random term that Frost includes, and that
omission must be justified and its effect on R re-derived. Re-run all of Section
5.5 and Tables/figure with corrected parameters, and reconcile against Frost
Table III as an external check.

**M3. Abstract (lines 29-33) and Conclusions (lines 38-44): the '30-45 percent'
reduction and 'simplest available closed-form' claims are not supported by the
manuscript's own results. [VERIFIED]**
The 30-45 percent figure is Frost's MIRIAD-specific headline, not a result the
present paper establishes. The manuscript's own run-in table (s2=1, t=1,
delta=0.5, J1=2) gives J0:0->2 reductions of 14, 28 and 76 percent across the
low/medium/high slope-variability columns; none of the three is in the 30-45
band, and the AD example (per M2) gives 99 percent. The 'simplest available
closed-form expressions' superlative is an unverifiable marketing claim and
conflicts with the user's own observation (Section 3.3.1) that the general case
'does not admit a simple closed-form expression'. *Remedy:* replace the abstract
numbers with figures actually produced by the corrected analysis, with the
parameter regime stated. Delete the 'simplest available' superlative.

**M4. Notation: silent reparameterisation and symbol collisions. [INSPECTED]**
(a) The abstract, model, and all prose use J0/J1/J2 for run-in/treatment/close,
but every code chunk and the R package use r/k/f. The reader cannot map the
implementation to the equations without guessing. (b) The symbol J is used both
for the total change-score dimension (J = J0 + J1 + J2, Eq. sigma context) and,
in Section 3.1, the per-participant scalar mathcal{J}_i; these are distinct but
visually close. (c) In Eq. (model) the letters a and b denote the random
intercept a_i and random slope b_i, then in Section 3.1 (line 360) a and b are
silently redefined as the variances sigma^2 and sigma_b^2 ('Let a = sigma^2 and
b = sigma_b^2'). This collides directly with the random effects a_i, b_i defined
20 lines earlier and matches the Mathematica scripts' convention rather than the
paper's. *Remedy:* unify J0/J1/J2 throughout including code; rename the
variance shorthands (the Mathematica c, d reparameterisation can stay in an
appendix); avoid reusing a, b.

**M5. Section 4 / model, lines 230-236: the common-close treatment-indicator
construction is asserted, not derived, and its claimed equivalence is not
checked. [INSPECTED, partially SUSPECTED]**
The model sets g_i = 0 during common close and states this is 'equivalent to
replacing the treatment indicator with g_i * I(1 <= j <= J_1)'. The resulting
gamma column (line 849-859, and design.R line 41) carries the value J_1 * t for
all common-close visits in the treatment arm, encoding the accumulated, frozen
treatment effect. This is a modelling choice with real consequences for what
gamma estimates (cumulative on-treatment divergence held constant off
treatment), and it presumes no carryover, no washout, and a perfectly
maintained level. None of these assumptions is stated. I verified the code
implements the stated column structure; I did not find a derivation showing the
construction yields an unbiased estimator of the intended estimand under the
generative model. *Remedy:* derive the common-close change-score mean explicitly
from the generative model (not by asserting g_i=0), state the carryover/washout
assumptions, and clarify the estimand. A delayed-start design (Liu-Seifert 2015;
Wang et al. 2019 delayed-start) makes the opposite assumption (effect may
attenuate); the contrast must be discussed.

**M6. The general case is claimed but never validated beyond the two closed-form
checkpoints. [INSPECTED]**
Section 3.3.1 (lines 792-806) states the general (J0, J1, J2) variance is
computed by the matrix routine and that J0=0 and J0=2 'serve as verification
checkpoints'. There is no independent validation of the general case, in
particular no check of the common-close path against a brute-force GLS
computation Var = [(X' Sigma^{-1} X)^{-1}]_{gamma,gamma} with Sigma assembled
directly (not via Woodbury), and no Monte Carlo confirmation. The Reproducibility
section (lines 1239-1243) mentions Monte Carlo benchmarking with a seed, but no
such benchmark appears in the manuscript. *Remedy:* add a direct
(non-Woodbury) GLS evaluation of Var(gamma) for several (J0, J1, J2) including
J2>0, and a small simulation, to confirm the Woodbury route and the design
matrices are correct in the common-close regime. The two existing checkpoints
both have J2=0 and so do not exercise the novel code path at all.

**M7. The X'R^{-1}X general formula (Eq. XRX_general, lines 760-773) is stated
without proof and one entry is asserted, not shown. [SUSPECTED, not fully
verified]**
The (delta, beta) cross term carries 2 D_{J0} D_{J1}/(J+1) with a sign argument
in prose (lines 778-780). I did not symbolically verify every entry of this 3x3
matrix; the numerical var_gamma_matrix agrees with the closed forms at the two
checkpoints, which exercises the assembled inverse but not each individual entry
of Eq. XRX_general in isolation. Because the manuscript presents Eq. XRX_general
as a result, each entry should be derivable. *Remedy:* either provide the
derivation (an appendix, or reference to the Mathematica output) or downgrade Eq.
XRX_general to 'the implementation computes', since it is not actually used in
closed form anywhere downstream.

**M8. Mathematica/MATLAB provenance does not match the manuscript's structure.
[VERIFIED]**
The CLAUDE.md and prior drafts reference frost311.m and frost312.m; the actual
files are woodbury_variance.m, blup_derivation.m, raw_outcome_model.m,
rcrm_raw_outcome.m. More importantly, woodbury_variance.m PART 2 derives the
one-run-in case (r=1, k=1, 2 change scores, R = a{{2,-1},{-1,2}}), which is
Frost's run-in design, not the manuscript's headline 'two run-in' (r=2, k=1)
case. So the symbolic script supporting the paper's Case 2 closed form is not
present; the r=2 closed form (Eq. var_r2k1) is asserted with an inline
cofactor/determinant sketch (lines 666-671) that I could not reproduce
symbolically from the script set. I did verify var_r2k1 numerically against the
matrix routine across 12 parameter combinations (max abs diff 3.6e-15), so the
formula is almost certainly correct; but the derivation trail is incomplete.
*Remedy:* supply the symbolic derivation for the r=2 case, or cite the exact
script that produces it. Fix the data-availability/reproducibility references
(the manuscript also calls the scripts 'MATLAB' in lines 1207, 1233, 1237 while
the syntax is Mathematica/Wolfram).

**M9. Verified-correct items (for the record). [VERIFIED]**
The following were re-derived or numerically reproduced and are correct: Eq.
var_r0k2 equals Frost Eq. 13 (V = 2 sigma_b^2 + sigma^2/t^2); the
Sherman-Morrison R^{-1} (Eq. Rinv); Z'R^{-1}Z for J=2 and J=3 (Eqs. ZRZ
context); var_gamma_matrix(1,1,0,.) equals Frost Eq. 16 exactly across 12
parameter sets; the Edland recovery (Eq. var_edland, Var = (2/n)[sigma_b^2 +
sigma^2/S_xx]) matches the matrix routine exactly for k in 2..4; var_r2k1
matches the matrix routine. The Woodbury application, per-participant additivity
(Eq. perpt_info), and the use of per-subject (not summed) outer products in the
W correction are correct and the W-construction subtlety (lines 636-655) is
handled properly in both the prose and design.R.

## 5. Minor comments

- Line 56 (knit hook) and lines 1245-1250: hard-coded absolute paths and a
  user-specific render hook will not survive submission; strip before sending.
- Line 274 (Eq. change_model) and line 185: the run-in index runs j = -J0,...,-1
  with randomization visit at j=0; but line 474 lists times '(-2t, -t, 0, t)'
  for J=3 change scores, listing four times for three change scores. The j=0
  baseline is differenced away and should not appear in the change-score time
  list. Tidy.
- Line 487-492: the sentence 'Z = t(1,1,1)' . (-2,-1,1)' . b_i' is dimensionally
  incoherent as written (a product of two column vectors and a scalar); the
  intended statement Z = t(-2,-1,1)' appears on the next line. Delete the
  malformed expression.
- Eq. H (lines 747-750): the displayed expression 'H = ab/(a + t^2 b * Z'R^{-1}Z
  * a/t^2)' is circular (it defines H in terms of Z'R^{-1}Z which itself contains
  a, t). Write H = b/eta_J using Eq. eta, which is already defined.
- citation style: the bib mixes two entry sets for the same works (e.g.
  frostOptimizingDesignClinical2008 and frost2008; wangTwoPeriodLinear2019 and
  wang2019; lairdlRandomEffectsModelsLongitudinal2023 and laird1982; iddi2022 is
  unique). De-duplicate; the second block (lines 463 onward) is an unprocessed
  Zotero import dump and includes ~50 entries never cited.
- lairdlRandomEffectsModelsLongitudinal2023 has corrupted author fields
  ('Lairdl, Nan M', 'Warel, James H') and a wrong year (the Laird-Ware paper is
  Biometrics 1982). The duplicate laird1982 entry is correct; cite that one.
- Tables 5.2-5.4: report which variance regime each column corresponds to in
  sigma_b/sigma terms in the caption (already partly done) and add the
  no-run-in baseline row's absolute n so reductions are auditable.
- Figure (relative efficiency, line 1065): axis labels use 'J_0' and 'J_2' as
  literal underscores in base-R plotting; they will not render as subscripts.
  Use expression() or switch to ggplot.
- Reproducibility (lines 1239-1243): claims a Monte Carlo benchmark with
  set.seed(20260418) that does not appear anywhere in the manuscript. Either add
  the benchmark or remove the claim (see M6).
- The DOI for frost2008 in the second bib block (line 483, number = {18}) gives
  the wrong issue; frostOptimizingDesignClinical2008 correctly has number = {19}.

## 6. Missing or questionable references

| Claim location | Issue | Suggested source (full citation) |
|---|---|---|
| Intro, run-in framework (lines 94-104) | Cites Frost 2008 but not its software realisation, which implements the same slope-power machinery and handles run-in | Nash S, Morgan KE, Frost C, Mulick A. Power and sample-size calculations for trials that compare slopes over time: introducing the slopepower command. Stata Journal. 2021;21(3):575-601. doi:10.1177/1536867X211045512. Plus the 2024 erratum, Stata Journal, doi:10.1177/1536867X241297951. |
| Intro / Section 3.3 baseline (lines 129-141) | longpower cited only in passing; should be the numerical benchmark for the boundary cases | Iddi S, Donohue MC. Power and sample size for longitudinal models in R: the longpower package and Shiny app. R Journal. 2022;14(1):264-282. (already in bib as iddi2022; elevate to a validation comparison) |
| Section 4 common close (lines 146-156) | Common-close framed without contrasting delayed-start, the natural comparator | Liu-Seifert H, Andersen SW, Lipkovich I, Holdridge KC, Siemers E. A novel approach to delayed-start analyses... PLoS One. 2015;10(3):e0119632. doi:10.1371/journal.pone.0119632. (in bib as liuseifert2015, uncited) |
| Section 6 J0=1 claim (lines 934-944) | Frison-Pocock cited for a result about a different estimator; either re-scope or remove | Frison L, Pocock SJ. Stat Med. 1997;16(24):2855-2872. (claim does not follow for the change-score GLS estimator; see M1) |
| 'Galbraith-Marschner' design guidance (if intended) | bib entry galbraith2002 lists Controlled Clinical Trials, which is correct; ensure it is actually cited, currently it is not | Galbraith S, Marschner IC. Guidelines for the design of clinical trials with longitudinal outcomes. Control Clin Trials. 2002;23(3):257-273. doi:10.1016/S0197-2456(02)00205-2. |
| 'lu2010' boundary-efficiency claim (line 120) | Year/journal need verification; bib has lu2010 as Biometrics 2010 'On efficiency of constrained longitudinal data analysis versus longitudinal ANCOVA' (this is real: Lu K, Biometrics 2010;66(3):891-896) | Confirmed correct as cited; the literature agent's flag was a false alarm. No change needed beyond confirming the DOI 10.1111/j.1541-0420.2009.01322.x. |
| Both hu2021 entries | The bib contains TWO distinct hu2021 keys (Biometrical Journal; Int J Biostatistics). BibTeX will silently drop one duplicate key. Both papers are real and distinct. | Rename keys (e.g. hu2021bimj and hu2021ijb) and cite each for its specific result. Hu N, Mackey HM, Thomas R. Biom J. 2021;63(4):806-824. AND Hu N, Mackey H, Thomas R. Int J Biostat. 2021;18(1):173-182. |
| Author self-positioning | Hu-Mackey-Thomas (2021) by the same author already gives closed-form random slope+intercept power; novelty delta not stated | (see M3, Significance) |

## 7. Suggestions for strengthening

1. Re-centre the novelty claim on the common-close phase. Derive a genuine
   closed form for at least J2=1 and J2=2 (the algebra is tractable for small
   J2 and is the contribution readers will care about), and validate it directly.
2. Add a direct-GLS and Monte Carlo validation harness covering J2>0, since no
   current check exercises the novel path (M6).
3. Provide an asymptotic statement: behaviour of Var(gamma) as J0, J1, J2 grow,
   and the limiting efficiency. Frost notes the run-in design's efficiency is
   unbounded in t (no asymptote), unlike the conventional design; an analogous
   characterisation for the common-close phase would add theoretical weight and
   help clear a higher-tier bar.
4. Map the model to Frost's three-component MIRIAD decomposition explicitly so
   the AD application is credible and reproduces Frost Table III as a sanity
   check (M2).
5. State the common-close estimand and its assumptions (no carryover, level
   maintenance) and contrast with delayed-start (M5).
6. Benchmark all boundary cases against longpower numerically and report the
   agreement, turning the prose claims of Section 3.3 into a table.
7. Unify J0/J1/J2 notation across prose, code, and package (M4).

## 8. Scope of review

**Verified (re-derived symbolically and/or reproduced numerically in R):**
Eq. var_r0k2 against Frost Eq. 13; var_gamma_matrix(1,1,0) against Frost Eq. 16
(12 parameter sets); Eq. var_r2k1 against the matrix routine (12 sets, max diff
3.6e-15); Eq. var_edland against the matrix routine (k=2..4); the R^{-1}
Sherman-Morrison form; Z'R^{-1}Z for J=2,3; the falsity of the J0=1 worsening
claim (exhaustive grid over k, sigma_b^2, t); the AD parameter
misappropriation against frost.pdf Table II; the abstract reduction figures
against the manuscript's own Tables; the self-contradictory AD sample sizes
(n=329 -> 2 -> 1). I reproduced the design-matrix construction (design.R) and the
Woodbury variance routine (variance.R) line by line.

**Inspected (read and judged, not independently re-derived):** Eq. XRX_general
entry by entry (M7); the common-close change-score derivation (M5); the
per-participant information form (Eq. perpt_info); the four Mathematica scripts
(woodbury_variance.m verified in part; blup_derivation.m, raw_outcome_model.m,
rcrm_raw_outcome.m read for consistency, the last confirms the RCRM benchmark
Var(beta1) = (sigma_e^2 + S v1)/(m S) reducing to Frost at S = 2t^2).

**Not checked:** full symbolic re-derivation of the r=2 cofactor/determinant
sketch (lines 666-671); the power.R/sample_size rounding behaviour at edge cases;
the AR1 extension (Paper 2) and allocation (Paper 3); the .csl/render pipeline;
any claim about cognitive-outcome variance regimes (line 1172) which is
qualitative.

**Literature searches run:** via a delegated web/PubMed agent, the following
were searched and the agent's findings incorporated: Frost-Kenward-Fox 2008;
Ard-Edland 2011; Liu-Liang 1997; Galbraith-Marschner 2002; Lu-Mehrotra-Liu 2009;
Lu 2010; Lu-Luo-Chen 2008; Hu-Mackey-Thomas 2021 (both journals, confirmed two
distinct real papers); Nash-Morgan-Frost-Mulick 2021 slopepower (and 2024
erratum); Wang et al. 2019 two-period run-in; delayed-start / randomized
withdrawal power literature; Iddi-Donohue 2022 longpower; and post-2015 run-in /
common-close developments. The Frost (2008) full text (docs/frost.pdf, 12 pages)
was read directly and used for the parameter and closed-form cross-checks.

---
*Rendered on 2026-06-13 at 18:18 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/analysis/report/01-runin-power/referee-report-2026-06-13.md*
