# Publication Review White Paper: Run-In Power Analysis Compendium
*Review date: 2026-08-16 16:04 PDT*

This is the third dated pub_review pass over the four manuscripts in
`analysis/report/`, following the whitepaper of 2026-08-07 (which
itself carried forward three prior referee reports dated 2026-06-13
through 2026-06-23). This pass re-reads all four manuscripts, checks
the repository's git history since the 2026-08-07 review session
(commits through 2026-08-15), and re-verifies the "still open" items
recorded in that whitepaper's revision history. Several quantitative
claims are re-executed against `R/*.R`; most findings in this round
rest on direct reading of the current `report.Rmd` files (inspected)
rather than code execution, and this is flagged per finding.
Epistemic labels: verified (code executed this session), inspected
(source read this session), inferred, unverified.

No repository changes were made during this review; only this
whitepaper was written.

---

## 1. Summary of the work under review

**Paper 1, `01-runin-power/report.Rmd` (1,527 lines).** Derives the
variance of the GLS treatment-effect estimator in a two-arm
longitudinal trial under a random-slopes model with change-score
parameterization, via the Woodbury identity, for three nested designs:
no run-in, two run-in observations, and a general
$(J_0, J_1, J_2)$ configuration with run-in and common close. The
MIRIAD example and the corrected residual variance are unchanged and
still verified. Since 2026-08-08 the manuscript has gained an explicit
paragraph distinguishing the common-close construction from the
delayed-start design of Liu-Seifert (2015), with the estimand
comparison explicitly and honestly deferred. The quantitative core
remains sound; the manuscript is materially unchanged in form since
the last review (still a two-layer drafting scaffold; see Major
Issue 1).

**Paper 2, `02-correlation/report.Rmd` (669 lines).** Compares AR(1)
and compound-symmetry (CS) residual correlation via a Woodbury/GLS
route and a moment-based route. The abstract and Discussion remain
correctly signed (AR(1) generally smaller than CS at matched nominal
$\rho$, anticonservative corner bounded near 6%), the two-route
equivalence claim remains correctly withdrawn, and the robustness
figure's clipping artifact (prior Minor Issue 5) is resolved. A titled
section, "Stratified variance for unequal visit patterns," still ends
in a bare formula shell and an explicit TODO comment with no content,
an issue this review treats as more serious than the prior "minor" TODO
classification because a numbered section with no delivered content is
a structural gap a referee will flag directly.

**Paper 3, `03-allocation/report.Rmd` (1,059 lines).** Poses the
fixed-budget allocation problem across run-in, treatment, and common
close phases. The Results-Discussion contradiction identified in the
2026-08-07 review (Major Issue 4) is resolved: the Discussion now
states the zero-run-in finding as the honest headline and explicitly
reconciles it with Paper 1's added-visit framing. The non-equal-spacing
section (Section 6) is now explicitly flagged as unreconciled with the
main change-score model rather than presented as complete, which is
an improvement in honesty though the underlying gap is unchanged.

**Paper 4, `04-gls-vs-ttest/report.Rmd` (1,224 lines).** Compares GLS,
ANCOVA, and averaged change-score estimators. The abstract now reports
the equal-estimand-gamma-scale rescaling from the 2026-08-08
remediation (ANCOVA RE 0.64-0.81 at $(2,2,0)$), which this review
re-executed and confirms to a first approximation (verified: RE 0.66-0.81
across signal-to-noise ratios 0.5-3 at $t=1$, matching the abstract's
range). The ordering proposition is stated and proved (inspected). The
title, however, still names only GLS and the averaged estimator,
omitting ANCOVA, the estimator the paper's own numbers favor as the
practical compromise; this specific finding was flagged as Minor Issue
7 on 2026-08-07 and remains unaddressed.

**Coherence across the set.** The shared computational engine, unified
bibliography (`references.bib`, symlinked from all four papers, 74
entries, no `zhao2022` duplicate), and standardized
`bookdown::pdf_document2` output format are all confirmed still in
place (inspected). The Paper 1 vs Paper 3 tension (run-in helps when
visits are added, never helps when visits are fungible) is now
explicitly reconciled in Paper 3's Discussion, closing the coherence
gap flagged in the prior review. The remaining coherence problem is
unchanged: all four abstracts and introductions still refer to "Paper
1 ... Paper 4 in this compendium" as internal cross-references rather
than citable units, and no paper's Data Availability statement gives a
repository URL or DOI.

## 2. Major issues

1. **All four manuscripts remain drafting scaffolds, not final prose.**
   (All papers, global.) `grep -c '{\.orig'` on the current files
   returns 49 / 17 / 33 / 23 blocks in Papers 1-4 respectively
   (verified by direct count), essentially unchanged from the
   2026-08-07 baseline of 48 / 17 / 31 / 23 (Paper 1 and 3 counts grew
   slightly, consistent with added content, not scaffold removal).
   Every body paragraph in every manuscript still renders inside a
   `::: {.orig data-latex=""}` block, which `claudecode.tex` still
   typesets as a gray, reduced-size quote environment (inspected,
   `claudecode.tex` lines 39-41). Since the last review a new artifact
   has appeared: each paper now has a companion `bullets.Rmd`
   supplement (424 / 177 / 312 / 226 lines) that duplicates the
   argument as boldface topic sentences plus bulleted itemize lists,
   explicitly labeled a "structured, section-by-section summary ...
   provided as a supplementary reading aid" (inspected,
   `01-runin-power/bullets.Rmd` line 38). This is a reasonable
   drafting aid but does not substitute for collapsing the main
   manuscript to final prose; no journal will referee a submission
   whose body text renders in a colored, reduced-size, explicitly
   provisional typeface. Remediation unchanged from 2026-08-07:
   promote the `.orig` text to plain paragraphs, delete the `orig`/
   `bullets`/`rgt` environments from `claudecode.tex`, and decide
   whether `bullets.Rmd` ships as supplementary material or is
   dropped.

2. **Paper 2's "Stratified variance for unequal visit patterns"
   section is an empty shell.** (`02-correlation/report.Rmd`, lines
   466-478.) The section states the general weighted-sum formula for
   heterogeneous per-participant $J_0$/$J_2$ counts and then, in place
   of a derivation, contains only `<!-- TODO: develop stratified
   computation for -->` / `<!-- both CS and AR(1) -->` (inspected)
   before jumping directly to numerical examples. This was listed as
   part of the general "(to be derived)" TODO inventory in the prior
   review's item 12(b); this review elevates it because the current
   text presents a section header and a formula as though the
   development were forthcoming within the document, which a referee
   will read as an unfulfilled promise rather than a scoping
   statement. Remediation: either derive the stratified variance for
   at least one structure (CS is a one-line extension of the existing
   moment-based route) or remove the section header and fold the
   formula into a Limitations paragraph that explicitly scopes the
   paper to balanced designs.

3. **Paper 4's title omits ANCOVA, the estimator the paper's own
   numbers recommend as the practical compromise.** (Title, page 1.)
   The current title is "GLS versus Averaged Change-Score Estimators
   for Clinical Trials with Run-In Observations" (inspected, line
   2-4). ANCOVA is one of the three estimators derived, appears in the
   abstract's headline result (RE 0.64-0.81), and is the estimator the
   2026-08-07 review's recommended framing designates as the
   practically relevant middle case; a reader selecting on title alone
   would not learn the paper addresses it. This exact gap was Minor
   Issue 7 on 2026-08-07 and is unresolved nine days later despite
   substantial other revision activity on this file (the abstract,
   ordering proposition, and all headline numbers were rewritten in
   the same window). Remediation: retitle to name all three
   estimators, e.g. "GLS, ANCOVA, and Averaged Change-Score Estimators
   for Clinical Trials with Run-In Observations."

4. **The common-close construction's carryover assumption is embedded
   in an equation but never stated as an assumption.**
   (`01-runin-power/report.Rmd`, lines 976-992.) The common-close
   change score is written with the treatment effect entering only
   through the accumulated term $\gamma g_i J_1 t$, i.e. the
   between-group offset acquired during treatment is carried forward
   unchanged (neither growing nor decaying) through the common-close
   phase. This is a substantive modeling choice, precisely the
   assumption the manuscript's own text (lines 176-187) identifies as
   central to distinguishing the design from delayed-start trials, yet
   the assumption itself is never named, justified, or flagged as a
   limitation; a search for "assum", "frozen", "carryover", and
   "level.maint" in the manuscript (inspected) finds the word
   "assumptions" only in the deferred delayed-start comparison
   paragraph and in two generic closing-Limitations sentences about
   equal spacing and no dropout (lines 1456-1463), neither of which
   names this construction. Remediation: add one sentence at the
   equation stating the assumption explicitly (constant post-treatment
   offset, no further divergence or convergence of group means) and
   cross-reference it from the delayed-start deferral paragraph.

5. **The r = 2 and general-J0 symbolic derivation trail is still
   absent.** (`analysis/scripts/woodbury_variance.m`.) The only
   Mathematica derivation file present in the working tree
   (`woodbury_variance.m`, 185 lines) covers exactly two cases: Part 1
   (no run-in, Frost Section 3.1.1) and Part 2 (one run-in
   observation, Frost Section 3.1.2); inspected in full, no r = 2 or
   general-$J_0$/$J_2$ block exists in this file, and a companion pair
   of files referenced by the repository's own `CLAUDE.md`
   (`frost311.m`, `frost312.m`, described as covering the two-run-in
   and generalized cases) is not present anywhere under
   `runinpower/analysis/`; the only surviving copy of `frost311.m` was
   found under `archive/frost311/`, outside the reviewed manuscript
   tree, and no `frost312.m` was found anywhere in the workspace. This
   status is unchanged since 2026-08-07 and the underlying gap
   (Major Issue 7 in that review, "no symbolic derivation trail")
   remains open; additionally, `CLAUDE.md`'s description of
   `analysis/data/raw_data/*.m` as the location of these files is now
   itself stale (the `.m` files live in `analysis/scripts/`, and
   `raw_data/` contains only a PDF), a repository-hygiene issue
   outside the manuscript scope but worth correcting alongside any
   remediation of this finding.

## 3. Minor issues

1. Paper 4's SUR (seemingly unrelated regression) subsection remains
   an explicitly flagged stub ("not developed in this paper ... is
   flagged as future work in the Discussion," lines 350-353,
   inspected) rather than being dropped as the 2026-08-07 review
   recommended. This is now honestly disclosed rather than a
   completeness overclaim, an improvement over the prior state, but
   the specification-without-derivation subsection still adds length
   without content; consider moving it to an appendix or dropping it
   at the next revision pass.
2. Paper 3's non-equal-spacing extension (Section 6, lines 480-518)
   is honestly flagged as unreconciled with the change-score model
   ("has not been established here," line 509) and its
   general-$(J_0,J_1,J_2)$ optimization-over-times subsection is a
   bare TODO comment (lines 514-518); both are consistent with the
   prior review's characterization but remain incomplete.
3. All four abstracts and introductions still cross-reference "Paper 1
   ... Paper 4 in this compendium" as internal text (inspected,
   e.g. `04-gls-vs-ttest/report.Rmd` line 25); unchanged from
   2026-08-07 Minor Issue 9.
4. No paper's Data Availability section gives a repository URL or DOI
   (checked via grep for `github.com`/`DOI`/`osf.io` in
   `01-runin-power/report.Rmd`; none found); the package `DESCRIPTION`
   also carries no `URL` or `BugReports` field. Unchanged from
   2026-08-07.
5. `CLAUDE.md`'s architecture description of the Mathematica
   derivation pipeline (`analysis/data/raw_data/*.m`, naming
   `frost311.m` and `frost312.m`) no longer matches the repository
   (files live in `analysis/scripts/`, and only `woodbury_variance.m`,
   `raw_outcome_model.m`, `blup_derivation.m`, and
   `rcrm_raw_outcome.m` are present; no `frost311.m`/`frost312.m`
   exist under `runinpower/`). This is project documentation, not a
   manuscript, but it will mislead a future contributor or reviewer
   trying to locate the r = 2 derivation flagged in Major Issue 5.
6. `share/MANIFEST.md` and the dated PDF snapshots under
   `analysis/report/share/` show three separate WIP render passes
   since the last review (2026-08-10, 2026-08-15 morning and
   afternoon) with no corresponding change in `.orig`/scaffold status;
   this indicates rendering and content-cleanup activity have been
   decoupled, worth noting only because it confirms the scaffold
   count in Major Issue 1 is not stale (the manuscripts have in fact
   been touched and re-rendered since 2026-08-08 without the scaffold
   being addressed).

## 4. What remains to be done

Ordered by importance for submission readiness. Items resolved since
2026-08-07 are omitted from this list; see Section 7 for the
resolution record.

**(a) Required for correctness**

1. Paper 2: derive or explicitly scope out the stratified-variance
   section rather than leaving a titled section with a TODO comment
   and no content (Major 2, new this round).
2. Paper 1: state the common-close carryover assumption explicitly at
   the point it is introduced (Major 4, new this round).
3. Paper 1: supply the r = 2 and general-entry symbolic derivation
   trail, or state plainly that the manuscript's general results rest
   on the numerical package routine (`var_gamma_matrix`) rather than a
   hand-derivable closed form, and adjust any language implying
   otherwise (Major 5, carried forward).

**(b) Required for acceptance**

4. Collapse the two-layer scaffold to final prose in all four papers
   and remove the `orig`/`rgt`/`bullets` environments from
   `claudecode.tex`; decide the disposition of the new `bullets.Rmd`
   companion documents (supplementary material or removed) (Major 1,
   carried forward, unresolved since 2026-06).
5. Paper 4: retitle to name ANCOVA (Major 3, carried forward from
   2026-08-07 Minor Issue 7, now elevated given it survived a full
   remediation pass on the same file).
6. Replace compendium-internal "Paper N" cross-references with
   citable preprint references, and add a repository URL/DOI to each
   Data Availability statement (Minor 3-4, carried forward).
7. Paper 3: complete or formally drop the non-equal-spacing extension
   (Section 6) and its optimization-over-times TODO (Minor 2, carried
   forward).
8. Paper 4: resolve or drop the SUR stub subsection (Minor 1, carried
   forward, disclosure improved but content still absent).

**(c) Desirable polish**

9. Correct `CLAUDE.md`'s stale description of the Mathematica
   derivation file locations and names (Minor 5, new this round).
10. Add a cross-paper consistency tinytest asserting the four papers'
    engines agree at shared parameter points (carried forward from
    2026-08-07, not re-checked this round; status unverified).

## 5. Recommended framing

The literature landscape and each paper's verified computational core
are unchanged since 2026-08-07, and the citation-engagement checks
performed this round (Paper 1's `liuseifert2015` delayed-start
paragraph, Paper 2's `lu2008`/`munoz1992` engagement, Paper 3's
optimal-design citations, Paper 4's `frison1992`/`frison1997`/
`liang1986` engagement, all inspected and confirmed substantive rather
than name-dropped) further support, rather than change, the prior
recommendations. This review carries forward the 2026-08-07
recommended framing for all four papers without modification, with one
addition specific to Paper 4.

**Paper 1.** Recommendation unchanged: frame around the common-close
phase as the headline result, with run-in as the established base
(Frost 2008, Wang 2019). Additional implication from this round:
because the r = 2 and general symbolic derivation trail remains
unsupplied (Major 5), the abstract and Methods should be explicit that
the general-$(J_0, J_2)$ result is established via the package's
numerical Woodbury routine and verified against direct GLS and Monte
Carlo, not via an independent hand derivation beyond r = 1; this is a
defensible and common practice in applied biostatistics papers but
should be stated rather than left to be inferred by a referee who goes
looking for the missing algebra.

**Paper 2.** Recommendation unchanged: retitle toward robustness of
run-in power formulas to serial correlation (framing ii from
2026-08-07), positioned against Lu-Luo-Chen (2008) and Frison-Pocock
(1992, 1997). The empty stratified-variance section (Major 2) should
either be completed as a short additional result under this framing
(a stratified CS variance is a direct extension of already-derived
moment formulas) or removed; under framing (ii) it is not essential to
the paper's core claim and can be cut without loss.

**Paper 3.** Recommendation unchanged: retitle around the fixed-budget
reconciliation with Paper 1 (framing ii from 2026-08-07). This
reconciliation is now actually written into the Discussion (Section
1, Major Issue 4 resolution), which strengthens the case for this
framing since the paper now demonstrably delivers the reconciliation
its title should promise.

**Paper 4.** Recommendation unchanged in substance (framing ii,
quantifying the efficiency cost of simpler estimators) but the title
implication from 2026-08-07 was not carried into this revision cycle
and should be treated as a required, not optional, title change at the
next pass: "GLS versus Averaged Change-Score Estimators" undersells a
paper whose principal quantitative contribution (RE 0.64-0.81 for
ANCOVA at $(2,2,0)$, re-verified this round) is about the
middle estimator, not the two named in the title.

## 6. Assessment

Compendium-level verdict: **major revision, not currently evaluable
for acceptance in any component; unchanged from 2026-08-07.**
Per-paper: Paper 1, major revision (quantitative core sound and
re-verified in part this round; the common-close estimand's central
assumption is still unstated and the symbolic derivation trail is
still incomplete); Paper 2, major revision, upgraded from the prior
"reject in present form" because both headline claims that motivated
that verdict are now correctly signed and the equivalence claim is
correctly withdrawn, but a titled section with no content is a new
finding a referee would flag at first read; Paper 3, major revision,
upgraded from the prior "reject in present form" because the
Results-Discussion contradiction that was the dispositive finding is
now resolved and the reconciliation with Paper 1 is written into the
text; Paper 4, major revision, upgraded from the prior "reject in
present form" because both abstract headlines are now correctly signed
and independently re-verified this round, but the title still does not
name the paper's central estimator nine days after every other part of
the same file was revised.

The trajectory since 2026-08-07 is real progress on exactly the
findings that whitepaper identified as dispositive: three of the four
"reject" verdicts move to "major revision" because the specific
false or self-contradicting claims that drove those verdicts
(Paper 2's inflation direction and equivalence claim, Paper 3's
Discussion-Results contradiction, Paper 4's both abstract headlines)
are now verifiably fixed rather than merely reworded. What has not
moved is the presentation layer common to all four papers: the
two-layer drafting scaffold (Major Issue 1) is unchanged in every
paper after two remediation sessions that touched the surrounding
prose extensively, and two specific, cheap, previously-identified
fixes (Paper 4's title, the compendium-internal cross-references) were
not applied despite substantial editing activity on the same files in
the same window. This pattern, content-level claims being reliably
fixed while structural and cosmetic items are reliably skipped even
when flagged repeatedly, is worth naming directly to the author: the
remaining path to submission readiness is now almost entirely
mechanical (title, scaffold removal, one missing sentence, one missing
paragraph, one missing citation location) rather than requiring new
derivation or new analysis, and the quantitative substance in all four
papers is, on the evidence re-checked this round, sound.

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
  (4,3,1)); disclosed in the manuscript. Resolved later the same day:
  all Paper 4 tables, figures, and the simulation were rescaled to
  the equal-estimand gamma scale (verified by end-to-end chunk
  execution); rescaled headlines are ANCOVA RE 0.64-0.81 at (2,2,0),
  MIRIAD N_ancova about 2.6x and N_averaged about 45x N_GLS at
  (4,3,1), and RE >= 0.95 in 12/115 cells in two characterized
  families. Still open (author-level): common-close estimand
  derivation and r = 2 symbolic trail (Paper 1), taxonomy content and
  stratified-variance TODO (Paper 2), Section 6 spacing model and
  optimal-design citations (Paper 3), cLDA/SUR/feasible-GLS
  simulation (Paper 4), scaffold collapse to final prose and
  compendium-internal cross-references (all).
- 2026-08-08 (continued): Literature-positioning pass, closing Major
  Issue 8. Verified by web search and added to the shared
  bibliography: tekle2008 (D-optimal cohort designs for linear
  mixed-effects models, Stat Med 27(14)), mentre1997 (optimal design
  in random-effects regression models, Biometrika 84(2)), munoz1992
  (damped-exponential correlation family, Biometrics 48(3)),
  atkinson1992 and fedorov1997 (optimum experimental design texts).
  Wired into the manuscripts: Paper 1 now cites nash2021
  (slopepower, with a stated scope difference), liuseifert2015 (a
  delayed-start contrast for the common-close phase, with the
  estimand comparison explicitly deferred), and galbraith2002; Paper
  2 cites munoz1992 as the family containing compound symmetry and
  AR(1) as boundary members; Paper 3 gained a paragraph positioning
  the allocation problem against the optimum-design literature and
  naming three features that distinguish it (constrained partition
  rather than free time points, a phase carrying its own nuisance
  parameter, and a c-optimality rather than D-optimality criterion);
  Paper 4 now cites frison1992 and frison1997 in the Introduction and
  at the Proposition with an explicit delimitation of what is new,
  liang1986 at the sandwich construction, and hu2021 in the run-in
  positioning. Also confirmed resolved: all four papers now symlink a
  single shared references.bib, closing the two-bibliography and
  two-key-scheme finding. Verified: every citation key in all four
  manuscripts resolves, no duplicate bib keys remain, and all four
  manuscripts execute end-to-end from purled chunks. Remaining open
  items are unchanged from the entry above, less the Frison-Pocock
  and optimal-design citation gaps now closed.
- 2026-08-16: Fresh review pass following one week of author activity
  (git log shows edits through 2026-08-15). Resolved since 2026-08-08:
  Paper 3's Results-Discussion contradiction (2026-08-07 Major Issue
  4) is fully reconciled in the Discussion text (verified by reading);
  Paper 2's robustness-figure clipping (Minor 5) is fixed (ylim now
  0-1.25, matching the caption's stated peak of 1.06); Paper 4's
  MVUE/BLUE terminology fix and title-vs-recommendation gap for the
  z/t critical values were confirmed already applied; the Paper 1 vs
  Paper 3 coherence tension is now explicitly reconciled in Paper 3's
  Discussion. Newly identified: Paper 2's "Stratified variance for
  unequal visit patterns" section is a titled section with only a
  formula and a TODO comment, no content (elevated to Major, not
  merely part of the general TODO inventory); Paper 1's common-close
  carryover assumption is embedded in an equation but never stated as
  an assumption (new Major finding); the Mathematica derivation file
  locations named in CLAUDE.md are stale relative to the actual
  repository layout (analysis/scripts/, not
  analysis/data/raw_data/), and frost311.m/frost312.m as named do not
  exist under runinpower/ (only an archived frost311.m survives
  outside the manuscript tree). Still open, unchanged since
  2026-08-08: the two-layer drafting scaffold across all four papers
  (`.orig` block counts 49/17/33/23, essentially unchanged from
  48/17/31/23, despite substantial intervening edits and the addition
  of new bullets.Rmd companion documents); the r = 2 / general-J0
  symbolic derivation trail (Paper 1); Paper 4's title still omitting
  ANCOVA (2026-08-07 Minor Issue 7, unresolved despite a full
  remediation pass touching the same file); compendium-internal
  "Paper N" cross-references and missing repository URL/DOI in all
  four Data Availability statements; Paper 3's non-equal-spacing
  Section 6 and Paper 4's SUR subsection, both now honestly disclosed
  as incomplete rather than misrepresented, but still incomplete.
  Re-verified this round by direct code execution: Paper 4's rescaled
  ANCOVA relative efficiency at (2,2,0) (RE 0.66-0.81 across
  signal-to-noise ratios 0.5-3, consistent with the abstract's stated
  0.64-0.81 range). Verdict change: Papers 2, 3, and 4 move from
  "reject in present form" to "major revision," reflecting that the
  specific false or self-contradicting claims that drove the prior
  reject verdicts are now verifiably corrected; Paper 1 remains major
  revision. Compendium-level verdict remains major revision overall.

---
*Rendered on 2026-08-16 at 16:04 PDT.*<br>
*Source: ~/prj/res/04-runin-power-analysis/runinpower/docs/pub_review_whitepaper_2026-08-16.md*
