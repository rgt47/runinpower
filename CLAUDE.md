# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`runinpower` is a research compendium (structured as an R package) for power analysis methodology in clinical trials with run-in observations. It focuses on optimal use of baseline rate information in longitudinal designs, particularly Alzheimer's disease trials with neuroimaging biomarkers. Author: Ronald G. Thomas (UCSD).

This is an analysis-focused project, not a distributable package. The R/ directory is intentionally minimal; the primary deliverable lives in `analysis/report/`.

## Development Commands

```bash
# Enter Docker container (recommended workflow)
make r

# Run tests
make test                # host: devtools::test()
make docker-test         # in container

# Full R CMD check
make docker-check

# Generate roxygen docs
make document

# Validate renv lockfile (default target, pure shell)
make check-renv

# Build Docker image
make docker-build

# Start RStudio Server (localhost:8787)
make rstudio
```

Single test file from R console: `testthat::test_file('tests/testthat/test-basic.R')`

## Architecture

**Two-tier reproducibility model:**

1. **Team layer (Dockerfile)** -- rocker/tidyverse:4.5.2 base with Quarto 1.6.43, system libs, Posit Package Manager for pre-compiled binaries
2. **Personal layer (renv.lock)** -- collaborative package accumulation; auto-snapshotted on container exit via `.Last()` hook in .Rprofile

**Container-first workflow:** Development happens inside Docker. Host R is optional. The Makefile and zzcollab CLI handle validation without requiring host R.

**Key directories:**

- `analysis/report/` -- RMarkdown report, references.bib, generated tables (tab1.tex/pdf)
- `analysis/data/raw_data/` -- Mathematica .m files (computational derivations)
- `analysis/data/derived_data/` -- processed outputs
- `R/` -- package functions (currently empty, pre-release)
- `tests/testthat/` -- testthat 3rd edition

## Environment Detection (.Rprofile)

The .Rprofile (v2.2.0) auto-detects context and adjusts behavior:

- **Container vs host:** `ZZCOLLAB_CONTAINER` env var; host skips renv, container uses Posit PM
- **CI guard:** `CI=true` disables auto-restore and auto-snapshot in GitHub Actions
- **IDE detection:** Suppresses noisy output in Neovim LSP and RStudio background sessions
- **Auto-snapshot on exit:** `.Last()` hook saves renv.lock changes when leaving container R

## CI/CD

GitHub Actions workflow (`.github/workflows/r-package.yml`) runs in rocker/tidyverse container: validates DESCRIPTION + renv.lock, caches renv packages, runs R CMD check and testthat.

## Dependencies

Minimal lockfile: renv 1.1.5 + testthat 3.3.1. Tidyverse stack comes from the Docker base image. Additional packages accumulate via renv as development proceeds.

## .Rbuildignore

The `analysis/`, `docs/`, `.github/`, Makefile, Dockerfile, and renv infrastructure are excluded from R package builds, separating the research compendium from any distributable package code.
