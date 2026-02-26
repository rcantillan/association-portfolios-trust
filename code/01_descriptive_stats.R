# ==============================================================================
# 01_descriptive_stats_independent.R — Descriptives from raw ELSOC longitudinal
# ------------------------------------------------------------------------------
# Independent of latent class / latent Markov steps.
# Inputs:
#   data/ELSOC_Long.RData (zip original loads: elsoc_long_2016_2022)
# Outputs (output/):
#   - summary_stats.csv
#   - membership_by_domain.csv
# ==============================================================================

source(here::here("code", "00_setup.R"))

# ---- Load raw data (robust to object name) ----
stop_if_missing(c(here::here("data", "ELSOC_Long.RData")))
obj_names <- load(here::here("data", "ELSOC_Long.RData"))

if ("elsoc_long_2016_2022" %in% obj_names) {
  raw <- get("elsoc_long_2016_2022")
} else if (length(obj_names) == 1) {
  raw <- get(obj_names[1])
} else {
  stop("Could not find elsoc_long_2016_2022. Objects: ", paste(obj_names, collapse = ", "))
}

dt_raw <- dplyr::as_tibble(raw)

# ---- Restrict to analysis waves used in paper (ola 1,3,6 -> 1,2,3) ----
dt <- dt_raw %>%
  dplyr::filter(ola %in% c(1, 3, 6)) %>%
  dplyr::mutate(
    id = idencuesta,
    ola = dplyr::case_when(ola == 1 ~ 1, ola == 3 ~ 2, ola == 6 ~ 3)
  )

# ---- Membership domains from c12_01..c12_08 ----
member_items <- paste0("c12_0", 1:8)
missing_member <- setdiff(member_items, names(dt))
if (length(missing_member) > 0) {
  stop("Missing membership items:\n- ", paste(missing_member, collapse = "\n- "))
}

# Original logic: if <2 => 0 else 1 (after missing cleanup)
dt <- dt %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(member_items), ~ to01(ifelse(to_na(.x) < 2, 0, 1)))) %>%
  dplyr::rename(
    nhg          = c12_01,
    religious    = c12_02,
    sport        = c12_03,
    charity      = c12_04,
    political    = c12_05,
    professional = c12_06,
    union        = c12_07,
    student      = c12_08
  )

domain_vars <- c("nhg","religious","political","union","professional","charity","sport","student")

# ---- Trust measures ----
# Generalized trust (zip original uses c02: 1->1, 2->0, 3->1)
if (!"c02" %in% names(dt)) stop("Missing c02 (generalized trust item).")

dt <- dt %>%
  dplyr::mutate(
    trust = dplyr::case_when(
      to_na(c02) == 1 ~ 1L,
      to_na(c02) == 2 ~ 0L,
      to_na(c02) == 3 ~ 1L,
      TRUE ~ NA_integer_
    )
  )

# Neighborhood trust: set in 00_setup.R (TRUST_NH_VAR), default "c03"
if (exists("TRUST_NH_VAR") && TRUST_NH_VAR %in% names(dt)) {
  dt <- dt %>% dplyr::mutate(trust_nh = to01(to_na(.data[[TRUST_NH_VAR]])))
} else {
  dt <- dt %>% dplyr::mutate(trust_nh = NA_integer_)
  warning("Neighborhood trust var not found or TRUST_NH_VAR not set. ",
          "Set TRUST_NH_VAR in 00_setup.R (e.g., 'c03' or 'c04').")
}

# ---- Covariates (minimal, consistent with your pipeline) ----
dt <- dt %>%
  dplyr::mutate(
    woman = dplyr::case_when(m0_sexo == 1 ~ 0L, m0_sexo == 2 ~ 1L, TRUE ~ NA_integer_),
    edad  = to_na(m0_edad),
    education = to_na(m01)
  )

# ---- Membership domains count (0–8) ----
dt <- dt %>%
  dplyr::mutate(membership_domains = rowSums(dplyr::across(dplyr::all_of(domain_vars), ~ to01(.x)), na.rm = TRUE))

# ---- Optional: balanced panel restriction (to match paper sample n=1,304) ----
# Toggle this to TRUE if you want descriptives for the balanced panel only.
USE_BALANCED_PANEL <- TRUE

if (USE_BALANCED_PANEL) {
  panel_n <- dt %>% dplyr::count(id)
  keep_ids <- panel_n %>% dplyr::filter(n == 3) %>% dplyr::pull(id)
  dt <- dt %>% dplyr::filter(id %in% keep_ids)
}

# ---- Education threshold (placeholder; adjust once you confirm m01 coding) ----
dt <- dt %>%
  dplyr::mutate(
    secondary_or_higher = dplyr::case_when(
      is.na(education) ~ NA_integer_,
      as.numeric(education) >= 3 ~ 1L,
      TRUE ~ 0L
    )
  )

# ----------------------------------------------------------------------
# (A) Main summary table by wave
# ----------------------------------------------------------------------
summary_tbl <- dt %>%
  dplyr::group_by(ola) %>%
  dplyr::summarise(
    N = dplyr::n(),
    age_mean = mean(as.numeric(edad), na.rm = TRUE),
    age_sd   = sd(as.numeric(edad), na.rm = TRUE),
    woman = mean(to01(woman), na.rm = TRUE),
    secondary_or_higher = mean(secondary_or_higher, na.rm = TRUE),
    generalized_trust   = mean(trust, na.rm = TRUE),
    neighborhood_trust  = mean(trust_nh, na.rm = TRUE),
    membership_domains_mean = mean(membership_domains, na.rm = TRUE),
    membership_domains_sd   = sd(membership_domains, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(summary_tbl, here::here("output", "summary_stats.csv"))
print(summary_tbl)

# ----------------------------------------------------------------------
# (B) Membership by domain (rates by wave)
# ----------------------------------------------------------------------
domain_tbl <- dt %>%
  dplyr::group_by(ola) %>%
  dplyr::summarise(
    N = dplyr::n(),
    dplyr::across(dplyr::all_of(domain_vars), ~ mean(to01(.x), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(domain_vars), names_to = "domain", values_to = "p_member") %>%
  dplyr::arrange(ola, domain)

readr::write_csv(domain_tbl, here::here("output", "membership_by_domain.csv"))
print(domain_tbl)

message("Balanced panel used: ", USE_BALANCED_PANEL)
if (exists("TRUST_NH_VAR")) message("Neighborhood trust item used: ", TRUST_NH_VAR)
message("NOTE: secondary_or_higher uses education >= 3; adjust to match m01 coding.")
