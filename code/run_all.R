# ==============================================================================
# run_all.R — One-command replication runner
# ------------------------------------------------------------------------------
# Run from project root:
#   source("code/run_all.R")
# ==============================================================================

source(here::here("code", "02_latent_markov.R"))
source(here::here("code", "01_descriptive_stats.R"))
source(here::here("code", "03_trust_models.R"))
source(here::here("code", "04_precarity_transitions.R"))

# Optional supplementary:
# source(here::here("code", "05_efa_cfa.R"))

message("Replication run completed. See /output for generated tables and figures.")
