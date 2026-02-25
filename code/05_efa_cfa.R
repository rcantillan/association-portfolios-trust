# ==============================================================================
# 05_efa_cfa.R — Optional EFA/CFA module (only if used in Supplementary)
# ------------------------------------------------------------------------------
# Input:
#   data/d1_wide.RData (expects object d1_wide)
# Output:
#   - output/efa_summary.txt
#   - output/cfa_summary.txt
# NOTE: This file is optional; it is not required for the main α/β/γ analyses.
# ==============================================================================

source(here::here("code", "00_setup.R"))

# Optional packages
req2 <- c("psych", "lavaan", "semTools")
installed <- rownames(installed.packages())
to_install <- setdiff(req2, installed)
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(req2, library, character.only = TRUE))

stop_if_missing(c(here::here("data", "d1_wide.RData")))
load(here::here("data", "d1_wide.RData"))  # loads d1_wide

d <- as_tibble(d1_wide)

# ---- Adapt item selection to your scale ----
trust_items <- d |> dplyr::select(dplyr::starts_with("trust_item"))

if (ncol(trust_items) < 3) {
  stop("Not enough trust_item* variables found in d1_wide. Update the selector in this script.")
}

# EFA
efa_fit <- psych::fa(trust_items, nfactors = 1, fm = "ml", rotate = "oblimin")
sink(here::here("output", "efa_summary.txt"))
print(efa_fit)
sink()

# CFA (update item names if needed)
item_names <- names(trust_items)[1:3]
model <- paste0("trust =~ ", paste(item_names, collapse = " + "))

cfa_fit <- lavaan::cfa(model, data = d, estimator = "WLSMV")
sink(here::here("output", "cfa_summary.txt"))
print(summary(cfa_fit, fit.measures = TRUE, standardized = TRUE))
print(semTools::reliability(cfa_fit))
sink()
