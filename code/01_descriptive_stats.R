# ==============================================================================
# 01_descriptive_stats.R — Descriptive tables used in the manuscript
# ------------------------------------------------------------------------------
# Inputs:
#   data/ELSOC_Long.RData
# Outputs (output/):
#   - summary_stats.csv
# ==============================================================================

source(here::here("code", "00_setup.R"))

stop_if_missing(c(here::here("data", "ELSOC_Long.RData")))
load(here::here("data", "ELSOC_Long.RData"))  # loads ELSOC_Long

dt <- as_tibble(ELSOC_Long)

# If your variable names differ, adapt these lines.
summary_tbl <- dt |>
  dplyr::group_by(ola) |>
  dplyr::summarise(
    N = dplyr::n(),
    generalized_trust = mean(to01(trust), na.rm = TRUE),
    neighborhood_trust = mean(to01(trust_nh), na.rm = TRUE),
    membership_domains = mean(memberships_count, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(summary_tbl, here::here("output", "summary_stats.csv"))
print(summary_tbl)
