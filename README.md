# Replication Package — Multiple Membership Configurations and Trust Formation 

This repository contains the replication materials for the manuscript:

**“Multiple membership configurations and trust formation: Structural precarity in highly unequal societies”**  
Roberto Cantillan, Gustavo Ahumada, Vicente Espinoza

## Overview

The paper examines whether the **configuration** of voluntary association memberships (portfolio positions) is associated with:
1) **generalized trust** and **neighborhood trust**, and
2) the **temporal stability** of cross-domain “bridging” portfolios (structural precarity).

Using three waves of ELSOC (2016, 2018, 2022; balanced panel), we:
- estimate a **latent Markov model** to classify respondents into three portfolio positions:
  - α = isolation
  - β = clustering (within-domain)
  - γ = bridging (cross-domain)
- analyze **transitions** across waves,
- estimate panel models (baseline RE probit) and report **Average Marginal Effects (AMEs)** for trust outcomes.

## Repository structure

- `code/`  
  Analysis scripts (R / Stata).
- `data/`  
  Input data objects used for replication (see Data availability).
- `output/`  
  Generated tables/figures used in the manuscript.
- `main.tex`, `preamble.tex`, `references.bib`  
  LaTeX manuscript source.

## Software requirements

- **R** (tested on R >= 4.1)
- Recommended: **RStudio**
- Key R packages: `LMest`, `tidyverse`, `panelr`, `survey`, `marginaleffects` (or equivalent AME tools)
- Optional: **Stata** (if using `code/estimation.do`)

> Reproducibility note: scripts are being standardized to run with **relative paths** from the project root
(using the `here` package). If you encounter hard-coded paths, see “Known issues” below.

## Data availability

This replication package uses the **Chilean Longitudinal Social Survey (ELSOC)** data.
ELSOC is administered by COES and may have redistribution constraints depending on the data agreement.

- If the data files are included in `data/`, replication should run out-of-the-box.
- If they are not included in your copy of the repository, users must obtain ELSOC from the official source
and place the files in `data/` following the file names indicated in the scripts.

## Reproducing the results

Run scripts in the following order:

### 1) Latent Markov model (portfolio positions α/β/γ)
- `code/002_long_latent_class.R`

Outputs (examples):
- latent position profiles
- transition matrix / transition plot(s)
- tables/figures saved to `output/`

### 2) Additional estimation / robustness (if applicable)
- `code/estimation.do` (Stata; optional)
- `code/01_efa_cfa.R` (if used for measurement checks)
- `code/descriptive_stats.Rmd` (descriptive tables)

### 3) Manuscript tables and figures
The main tables/figures are stored in `output/` and are referenced in the LaTeX manuscript:
- `output/plot_transition.png`
- `output/emd_plot.png`
- `output/m3_fit_table.tex`, etc.

## Expected outputs

After running the scripts, the repository should contain in `output/`:
- Model fit table for the latent Markov solution (BIC/AIC/logLik)
- Position profiles (α/β/γ)
- Transition probabilities (and plot)
- Trust models with AMEs for generalized and neighborhood trust
- Supplementary tables (e.g., multinomial predictors of position membership)

## Citation

If you use these materials, please cite the manuscript and this repository.  
A `CITATION.cff` file will be provided upon journal submission.

## Contact

Roberto Cantillan — roberto.cantillan@...
