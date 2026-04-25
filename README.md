# Replication Package — Why Associational Life Fails to Build Generalized Trust under Inequality

This repository contains the replication materials for the manuscript:

**"Structural Precarity and the Civil Society–Trust Gap: Portfolio Configurations and Durability Under Inequality"**  

---

## Summary of the contribution

The paper argues that the relationship between voluntary associations and generalized trust is not well captured by additive participation measures. What matters is the **configuration of membership portfolios across organizational domains** and their **durability over time** under conditions of material inequality.

Using three waves of ELSOC (2016, 2018, 2022; balanced panel, *n* = 1,304), we:

1. Estimate a **latent Markov model** to classify respondents into three portfolio positions:
   - **α (isolation)** — minimal organizational embeddedness
   - **β (clustering)** — memberships concentrated within a single domain
   - **γ (bridging)** — memberships spanning multiple heterogeneous domains
2. Analyze **position transitions** to test structural precarity: γ is the least stable position (*P*(stay) = 0.58), and exits from γ are socioeconomically patterned
3. Estimate panel random-effects probit models for **generalized trust** and **neighborhood trust**, reporting **Average Marginal Effects (AMEs)**
4. Document the **double reversal**: γ > α > β for generalized trust; γ > β > α for neighborhood trust — the key theoretical prediction distinguishing a configurational account from additive alternatives

---

## Repository structure

```
code/
├── 00_setup.R                        # Shared config, STATE_MAP, helper functions, package loading
├── run_all.R                         # Single entry point — runs full pipeline in order
├── main/
│   ├── 01_descriptive_stats.R        # Sample description, summary statistics
│   ├── 02_latent_markov.R            # LMM estimation (K=1..5), profile extraction, state mapping
│   ├── 03_trust_models.R             # RE probit trust models, AMEs, γ−β contrast
│   ├── 06_figures_plots.R            # All main manuscript figures
│   └── estimation.do                 # Stata cross-check (optional)
└── SI/
    ├── 03a_trust_sensitivity.R       # S1–S5: CRE/Mundlak, Oster δ, cross-lagged,
    │                                 #         entropy-weighted, soft assignment
    ├── 03b_precarity_transitions_SES.R  # H3: exit rates, SES predictors, bootstrap (B=999)
    ├── 03c_active_only_sensitivity.R    # S6: active-only membership sensitivity
    ├── 03_H4_connectivity_best.R        # H4: portfolio connectivity (Jaccard)
    ├── 04_SI_tables.R                   # All SI LaTeX tables
    └── 07_replicate_stata_models.R      # Stata/R cross-validation

data/                                 # Input data and intermediate objects (see below)
output/                               # Generated tables (.txt) and figures (.pdf) — main + SI
```

> **Paths:** All R scripts use relative paths via the `{here}` package. Run from the repository root.  
> **Membership coding:** `MEMBER_CODE_LOGIC = "any_member"` (c12 ≥ 2) is the main analysis specification. Active-only sensitivity (c12 == 3) is in `code/SI/03c_active_only_sensitivity.R`.  
> **State mapping (canonical):** α = State 3 (Isolation), β = State 1 (Clustering), γ = State 2 (Bridging). This mapping is defined once in `00_setup.R` and used throughout.

---

## Requirements

### R environment

- **R ≥ 4.1** recommended
- Key packages (installed automatically if missing by `00_setup.R`):
  - `LMest` — latent Markov model estimation (EM algorithm)
  - `lme4` — random-effects probit (`glmer`, `nAGQ = 12`)
  - `marginaleffects` — average marginal effects
  - `sandwich`, `clubSandwich` — cluster-robust standard errors
  - `data.table`, `tidyverse`, `here`, `broom`, `ggplot2`, `patchwork`

### Optional: Stata

A Stata do-file is included as a cross-check of the RE probit + margins results:
- `code/main/estimation.do`
- Note: `xtprobit` uses Gauss-Hermite quadrature (nAGQ = 9) by default; R's `glmer` defaults to Laplace (nAGQ = 1). Set `nAGQ = 12` in R to match Stata. Coefficient differences of ~8× between defaults are expected and documented in the cross-validation script.

---

## Data availability (ELSOC)

This project uses **ELSOC (Chilean Longitudinal Social Survey)** data, administered by the Center for Social Conflict and Cohesion Studies (COES). Redistribution is subject to COES data use terms.

- If the required data files are present in `data/`, replication runs out-of-the-box via `run_all.R`.
- If not included, obtain ELSOC from the official source and place the required object in `data/`.

**Official data source:** [https://coes.cl/encuesta-panel/](https://coes.cl/encuesta-panel/)

### Expected input

The pipeline expects:

- `data/ELSOC_Long.RData` — containing an object named `ELSOC_Long`

Required columns: `id`, `ola` (wave: 1 = 2016, 2 = 2018, 3 = 2022), the 8 membership domains (`nhg`, `religious`, `political`, `union`, `professional`, `charity`, `sport`, `student`), trust outcomes (`trust`, `trust_nh`), and covariates (`edad`, `woman`, `education`, `employment`).

> **Missing value convention:** ELSOC codes item refusals as −888 and don't-know responses as −999. The pipeline recodes these explicitly to `NA` before analysis. See `00_setup.R` for details.

---

## How to reproduce all results

From the repository root, run:

```r
source("code/run_all.R")
```

This executes the full pipeline in order:
1. `00_setup.R` — environment setup
2. `main/01` → `main/03` → `main/06` — main manuscript results and figures
3. `SI/03a` → `SI/03b` → `SI/03c` → `SI/04` — supplementary analyses and tables

All outputs (tables as `.txt`, figures as `.pdf`) are written to `output/`.

To reproduce a single stage, source the relevant script directly — each script loads its dependencies via `here::here()` and the shared config from `00_setup.R`.

---

## Key numerical results (for verification)

| Result | Value |
|---|---|
| Balanced panel N | 1,304 individuals; 3,891 person-waves |
| γ (Bridging) prevalence | ~11% of person-wave observations |
| γ persistence probability (model-implied) | 0.31 |
| β persistence probability (model-implied) | 0.61 |
| AME of high education on P(exit from γ) | −0.220 (SE = 0.065, p < .001) |
| AME β vs. γ — generalized trust (main model) | −0.056 (SE = 0.018, p < .01) |
| AME α vs. γ — generalized trust (main model) | −0.033 (SE = 0.017, p < .10) |
| AME α vs. γ — neighborhood trust (main model) | −0.080 (SE = 0.033, p < .05) |
| CRE/Mundlak AME β vs. γ — generalized trust | −0.081 (SE = 0.023, p < .001) |

---
