# ==============================================================================
# run_all.R — Master script (one-command replication)
# ==============================================================================

source(here::here("code", "00_setup.R"))

message("1) Descriptives...")
source(here::here("code", "01_descriptive_stats.R"))

message("2) Latent Markov model...")
source(here::here("code", "02_latent_markov.R"))

message("3) Trust models + AMEs + γ−β contrast...")
source(here::here("code", "03_trust_models.R"))

message("4) Precarity transitions (exit γ)...")
source(here::here("code", "04_precarity_transitions.R"))

message("DONE. Check /output for tables and figures.")