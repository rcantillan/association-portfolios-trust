# Replication Package — Multiple Membership Configurations and Trust Formation (ELSOC, Chile)

This repository contains the replication materials for the manuscript:

**“Multiple membership configurations and trust formation: Structural precarity in highly unequal societies”**  
Roberto Cantillan, Gustavo Ahumada, Vicente Espinoza  

## Summary of the contribution

The paper argues that the relationship between voluntary associations and trust is not well captured by additive participation measures. What matters is the **configuration of membership portfolios across organizational domains** and their **durability over time**.

Using three waves of ELSOC (2016, 2018, 2022; balanced panel), we:
1) estimate a **latent Markov model** to classify respondents into three portfolio positions  
   - **α (isolation)**, **β (clustering)**, **γ (bridging)**  
2) analyze **position transitions** to test structural precarity (lower stability of γ)  
3) estimate panel models for **generalized trust** and **neighborhood trust**, reporting **Average Marginal Effects (AMEs)**  
4) report the key inferential contrast for SSR: **γ − β** on generalized trust (Wald / diff-in-AME)

---

## Repository structure

- `code/`  
  Main analysis scripts (R). Includes a single entry point `run_all.R`.
- `data/`  
  Input data objects required to run the pipeline (see “Data availability”).
- `output/`  
  All generated tables and figures used in the manuscript and Supplementary.
- `manuscript/` (optional if present)  
  LaTeX sources (`main.tex`, `preamble.tex`, references).

> All R scripts use **relative paths** via the `{here}` package. Run scripts from the repository root.

---

## Requirements

### R environment
- **R >= 4.1** recommended  
- Key packages used (installed automatically if missing by `code/00_setup.R`):
  - `data.table`, `tidyverse`, `here`, `LMest`, `marginaleffects`, `lme4`, `broom`, `ggplot2`

### Optional: Stata
A Stata do-file is included as an optional cross-check of the RE probit + margins results:
- `code/estimation.do`

---

## Data availability (ELSOC)

This project uses **ELSOC (Chilean Longitudinal Social Survey)** data administered by COES. Redistribution may be restricted depending on the usage agreement.

- If the required data files are included in `data/`, replication runs out-of-the-box.
- If not included, users must obtain ELSOC from the official source and place the required objects in `data/`.

### Expected input files
The pipeline expects at minimum one of the following (depending on your setup):

- `data/ELSOC_Long.RData` containing an object named `ELSOC_Long`

The object is expected to include (at least) the following columns:
`id`, `ola`, the 8 membership domains (`nhg`, `religious`, `political`, `union`, `professional`, `charity`, `sport`, `student`), trust outcomes (`trust`, `trust_nh`), and covariates (`edad`, `woman`, `education`, `employment`).

---

## How to reproduce all results (recommended)

From the repository root, run:

```r
source("code/run_all.R")