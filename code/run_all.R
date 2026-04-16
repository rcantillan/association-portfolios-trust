# ==============================================================================
# run_all.R — Master script (one-command replication)
# ==============================================================================
# Structure:
#   code/00_setup.R       — shared configuration (MEMBER_CODE_LOGIC, STATE_MAP, etc.)
#   code/main/            — main manuscript analysis
#   code/SI/              — supplementary information analysis
#
# Usage:
#   Rscript code/run_all.R          # full replication
#   source("code/run_all.R")        # from R console at project root
# ==============================================================================

source(here::here("code", "00_setup.R"))

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

message("\n[1/4] Descriptive statistics...")
source(here::here("code", "main", "01_descriptive_stats.R"))

message("\n[2/4] Latent Markov model (K=3, any_member)...")
source(here::here("code", "main", "02_latent_markov.R"))

message("\n[3/4] Trust models + AMEs + gamma-beta contrast...")
source(here::here("code", "main", "03_trust_models.R"))

message("\n[4/4] Main figures...")
source(here::here("code", "main", "06_figures_plots.R"))

# ==============================================================================
# SUPPLEMENTARY INFORMATION
# ==============================================================================

message("\n[SI-1] H3: Structural precarity / exit from gamma...")
source(here::here("code", "SI", "03b_precarity_transitions_SES.R"))

message("\n[SI-2] H4: Portfolio connectivity (Jaccard)...")
source(here::here("code", "SI", "03_H4_connectivity_best.R"))

message("\n[SI-3] Sensitivity analyses S1-S5 (trust models)...")
source(here::here("code", "SI", "03a_trust_sensitivity.R"))

message("\n[SI-4] S7: Active-only LMM sensitivity...")
source(here::here("code", "SI", "03c_active_only_sensitivity.R"))

message("\n[SI-5] SI tables (LaTeX)...")
source(here::here("code", "SI", "04_SI_tables.R"))

message("\n[SI-6] Stata cross-check replication...")
source(here::here("code", "SI", "07_replicate_stata_models.R"))

message("\nDONE. Check output/ for all tables and figures.")
