# pub_review Remediation Log
*2026-08-16 16:23 PDT*

Remediation pass against
`docs/pub_review_whitepaper_2026-08-16.md`, addressing the
"What remains to be done" checklist (Section 4) and the Major
Issues (Section 2). This log records what was actually changed;
it does not replace the whitepaper, which remains the review
record and was not edited.

All four manuscripts were re-purled and executed end-to-end
after these edits (`knitr::purl()` + `source()` on each
`report.Rmd`), and the `tinytest` suite (68 assertions across 8
files) was re-run; both pass. See "Tests run" below for exact
commands and output.

## Fixed

1. **Major 4 / checklist (a)2 -- Paper 1's common-close carryover
   assumption was never stated as an assumption.**
   `analysis/report/01-runin-power/report.Rmd` (lines ~988-1010,
   ~186-189). Added a "Carryover assumption" paragraph at the
   point the common close change score is defined, naming the
   assumption explicitly (frozen, non-decaying between-group
   offset through the common close phase) and its consequence
   for the estimand if post-cessation trajectories in fact
   converge or diverge. Cross-referenced from the delayed-start
   deferral paragraph in the Introduction. `[applied, unverified
   for prose correctness; verified that the manuscript still
   purls and executes end-to-end after the edit]`.

2. **Major 5 / checklist (a)3 -- Paper 1's r=2 and general-J0
   symbolic derivation trail is absent, but the manuscript's
   Reproducibility section implied full symbolic coverage.**
   `analysis/report/01-runin-power/report.Rmd` (Reproducibility
   section, ~line 1519; Limitations paragraph, ~line 1484).
   Rewrote the Reproducibility paragraph to state precisely
   which cases `woodbury_variance.m` covers ($J_0 \in \{0,1\}$
   only) and that the general and $J_0=2$ results come from the
   package's numerical Woodbury routine
   (`var_gamma_matrix()` in `R/variance.R`), validated against a
   direct GLS computation and a Monte Carlo benchmark, not an
   independent symbolic derivation. Added a matching sentence to
   the Limitations paragraph. `[verified: grepped R/variance.R
   to confirm var_gamma_matrix() location; confirmed
   woodbury_variance.m contains only Part 1 (J0=0) and Part 2
   (J0=1) by prior whitepaper inspection, re-confirmed by file
   read]`.

3. **Major 2 / checklist (a)1 -- Paper 2's "Stratified variance
   for unequal visit patterns" section was a titled section with
   only a formula and a TODO comment.**
   `analysis/report/02-correlation/report.Rmd` (was lines
   466-478). Removed the empty section header and formula-only
   shell; folded the content into a new "Scope: balanced
   designs only" paragraph appended to the Discussion, stating
   the weighted-sum stratification idea, noting it is a
   straightforward extension for the compound-symmetry route but
   has not been derived or verified for either structure, and
   that unbalanced-visit designs are explicitly out of scope
   rather than promised. `[applied, verified: report.Rmd purls
   and executes end-to-end after the edit]`.

4. **Major 3 / checklist (b)5 -- Paper 4's title omitted ANCOVA.**
   `analysis/report/04-gls-vs-ttest/report.Rmd` (title, line
   2-3) and `analysis/report/04-gls-vs-ttest/bullets.Rmd` (title
   and body cross-reference to the main title, lines 2 and 41).
   Retitled to "GLS, ANCOVA, and Averaged Change-Score Estimators
   for Clinical Trials with Run-In Observations" in both files.
   The corresponding `.tex` files (`report.tex`, `bullets.tex`)
   are render artifacts and will pick up the new title on the
   next `bash tools/render.sh` pass; they were not hand-edited.
   `[verified: grep confirms no remaining "GLS versus Averaged"
   occurrences in the .Rmd sources]`.

5. **Minor 4 / checklist (b)6, partial -- no repository URL/DOI
   in Data Availability statements or `DESCRIPTION`.**
   `DESCRIPTION` (added `URL`/`BugReports`),
   `analysis/report/{01,02,03,04}-*/report.Rmd` (Data
   availability statement in each). Confirmed a real GitHub
   remote exists (`git remote -v` ->
   `git@github.com:rgt47/runinpower.git`) and that
   `https://github.com/rgt47/runinpower` returns HTTP 200
   (public), then added the URL to `DESCRIPTION` and to each
   paper's Data Availability statement. `[verified: curl
   confirmed HTTP 200; grep confirms the URL now appears in all
   four report.Rmd files and DESCRIPTION]`.

6. **Minor 2 / checklist (b)7 -- Paper 3's non-equal-spacing
   Section 6 optimization-over-times subsection was a bare TODO
   comment.** `analysis/report/03-allocation/report.Rmd` (was
   lines 512-518). Removed the "## Optimization over observation
   times within phases" subsection heading (which promised
   content it did not deliver) and replaced the bare `<!-- TODO
   -->` comment with an explicit "Scope" paragraph stating the
   general free-time-point optimization is not addressed and is
   left as future work, with an informal pointer to the
   optimal-time-point literature. The "General time points"
   subsection with its actual (raw-outcome-model) derivation was
   left in place unchanged, since it has real content, not a
   TODO. `[applied, verified: report.Rmd purls and executes
   end-to-end after the edit]`.

7. **Minor 5 / checklist (c)9 -- `CLAUDE.md`'s description of the
   Mathematica derivation pipeline was stale.** `CLAUDE.md`
   (Architecture section, "Derivation pipeline" paragraph).
   Corrected the path (`analysis/scripts/`, not
   `analysis/data/raw_data/`), replaced the nonexistent
   `frost311.m`/`frost312.m` file names with the actual present
   file (`woodbury_variance.m`) and its true scope ($J_0 \in
   \{0,1\}$ only), and noted where the archived `frost311.m`
   survives and that no `frost312.m` exists anywhere in the
   workspace. `[verified: file read confirms the corrected paths
   and file names match the actual repository layout]`.

8. **Checklist (c)10 -- no cross-paper consistency test existed.**
   New file
   `inst/tinytest/test_cross_paper_consistency.R` (7 new
   assertions). Added tests that (a) Paper 2's AR(1) engine
   (`var_gamma_ar1`) at rho=0 matches Paper 1's Woodbury engine
   (`var_gamma_matrix`) exactly, across four
   $(J_0,J_1,J_2)$ configurations, since the AR(1) residual
   covariance algebraically reduces to compound symmetry at
   rho=0; (b) Paper 3's `marginal_reduction()` reproduces the
   corresponding `var_gamma_matrix()` difference; (c) Paper 4's
   `re_ancova_vs_gls()` reproduces the corresponding
   `var_gamma_matrix()`/`var_gamma_ancova()` ratio and never
   exceeds 1 (GLS at least as efficient as ANCOVA). `[verified:
   ran `tinytest::run_test_dir("inst/tinytest")` after
   `pkgload::load_all(".")`; all 68 assertions across 8 test
   files pass, including these 7 new ones]`.

## Deferred

- **Major 1 / checklist (b)4 -- collapse the two-layer `.orig`/
  `bullets`/`rgt` scaffold to final prose in all four manuscripts
  and remove the corresponding environments from
  `claudecode.tex`.** Out of budget: this requires rewriting
  approximately 122 `.orig` blocks (49+17+33+23) across
  ~4,500 combined manuscript lines from a colored-quote
  provisional format to final plain prose, a task that is
  large, high-risk to do quickly (re-flowing prose without
  altering technical content), and explicitly time-boxed out of
  this remediation pass per the task's budget instructions. No
  partial or cosmetic fix was attempted, since a half-collapsed
  scaffold (some sections converted, others not) would be worse
  than a consistently-scaffolded document for a referee to read.
  Deferred to the author. No single command reproduces this; it
  requires manual editorial work file-by-file plus removing the
  `orig`/`rgt`/`bullets` environment definitions from each
  `claudecode.tex` (four files) and a decision on whether
  `bullets.Rmd` ships as supplementary material.
- **Minor 3 / checklist (b)6, remaining half -- replace
  compendium-internal "Paper N" cross-references with citable
  preprint references.** Not fixable without a decision only the
  author can make: no preprint of any of the four papers has
  been posted (no arXiv/medRxiv/bioRxiv identifier exists to
  cite), so replacing "Paper 1...Paper 4 in this compendium"
  with `@citekey` references would require either (a) posting
  preprints first and obtaining DOIs, or (b) the author deciding
  on a different citation mechanism (e.g., a working-paper
  series identifier). Fabricating placeholder DOIs was avoided
  per the hard constraint against fabricated results. The
  repository-URL half of this item (Minor 4) was completed (see
  Fixed #5).
- **Minor 1 / checklist (b)8 -- Paper 4's SUR stub subsection.**
  Reviewed (`analysis/report/04-gls-vs-ttest/report.Rmd`, "##
  SUR" subsection, lines ~331-353, and the correctly
  cross-referenced Discussion paragraph, "Section 2.4"). Left
  unchanged: the subsection is already honestly disclosed as
  unfinished ("not developed in this paper... flagged as future
  work"), and the Discussion's cross-reference to "Section 2.4"
  is in fact correct given the manuscript's actual heading
  numbering (verified by counting `#`/`##` headings), so this is
  not a correctness defect, only an acceptance-tier length/
  content question. Completing the SUR variance derivation under
  time pressure risked introducing an unverified or incorrect
  formula, which the task's constraints explicitly prohibit;
  dropping it outright was judged lower-value than the
  correctness-tier and cheap acceptance-tier fixes actually
  completed. Recommendation for the author: either derive the
  SUR variance (a modest but nontrivial GLS-on-stacked-equations
  calculation) or delete the "## SUR" subsection and its
  Discussion sentence at the next revision pass.
- **Full-scale simulation reruns.** Not needed this pass: the
  whitepaper's dead/empty-simulation-output defect class did not
  apply this round (no NA estimates or crashed runs were
  flagged in the 2026-08-16 whitepaper), so no simulation rerun
  was required. All four manuscripts were confirmed to execute
  end-to-end from purled chunks with no errors (see Tests run).
- **Testing-framework mismatch in `CLAUDE.md`** (new issue, not
  a whitepaper item; see below) was left uncorrected beyond the
  Derivation-pipeline paragraph already fixed, to avoid scope
  creep beyond the whitepaper's checklist within the time
  budget.

## New issues found while fixing

- `CLAUDE.md`'s "Development Commands" and "Single test file"
  sections describe `testthat` (`devtools::test()`,
  `testthat::test_file('tests/testthat/test-basic.R')`), but the
  actual test suite lives entirely under `inst/tinytest/` and
  uses `tinytest`; there is no `tests/testthat/` directory in
  this repository. This will mislead a contributor trying to run
  or add tests. Not fixed in this pass (out of the whitepaper's
  checklist scope); flagged for a future documentation pass.
- The `inst/tinytest/` file naming (`test_allocation.R`,
  `test_ancova.R`, etc., underscore-separated) does not match
  this user's stated global convention of `test-*.R`
  (hyphen-separated, per `~/.claude/CLAUDE.md`). This predates
  the current remediation and was not renamed, since renaming
  test files is unrelated to any whitepaper finding and outside
  this task's scope; flagged for author awareness.
- Several files were already modified in the working tree before
  this remediation session began (`analysis/report/*/bullets.pdf`,
  `bullets.tex`, and three of the four `bullets.Rmd` files, plus
  `analysis/report/share/MANIFEST.md` and new WIP PDFs under
  `analysis/report/share/`), consistent with the whitepaper's
  Minor Issue 6 note about decoupled render passes. These
  pre-existing changes were left untouched except for the
  title-string edit made to `04-gls-vs-ttest/bullets.Rmd` in Fix
  #4 above.

## Tests run

```
Rscript -e 'pkgload::load_all("."); tinytest::run_test_dir("inst/tinytest")'
# All ok, 68 results (0.8s)
```

End-to-end execution check (purl + source, all four manuscripts):

```r
pkgload::load_all(".", quiet = TRUE)
for (f in c("analysis/report/01-runin-power/report.Rmd",
            "analysis/report/02-correlation/report.Rmd",
            "analysis/report/03-allocation/report.Rmd",
            "analysis/report/04-gls-vs-ttest/report.Rmd")) {
  out <- tempfile(fileext = ".R")
  knitr::purl(f, output = out, documentation = 0, quiet = TRUE)
  source(out, local = new.env())
}
# All four: EXEC OK, no errors
```

`.Rmd` syntax check (`knitr::purl()` + `parse()`, no execution):
all four files parse without error.

No `.Rmd` was rendered to PDF through `tools/render.sh` in this
session (rendering was optional per the task instructions and
was skipped to stay within budget); a render is still recommended
before the next external submission or referee pass to confirm
the LaTeX output (in particular the retitled Paper 4 and the
edited Paper 1/Paper 2/Paper 3 sections) typesets cleanly. Exact
commands for the user:

```bash
bash tools/render.sh analysis/report/01-runin-power/report.Rmd
bash tools/render.sh analysis/report/02-correlation/report.Rmd
bash tools/render.sh analysis/report/03-allocation/report.Rmd
bash tools/render.sh analysis/report/04-gls-vs-ttest/report.Rmd
```
