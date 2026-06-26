# ── 01_descriptive_stats.R — Descriptives and balanced panel construction ────
#
# Balancing replicates Stata estimation.do:
#   (a) complete.cases() on all model covariates (equiv. probit e(sample))
#   (b) keep only ids observed in all 3 waves
#
# Outputs:
#   output/summary_stats_balanced.csv
#   output/summary_stats_full.csv
#   output/membership_by_domain.csv
#   data/dt_analysis.rds

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(data.table)
})

# ── 1. Load data ──────────────────────────────────────────────
stop_if_missing(c(here::here("data", "ELSOC_Long.RData")))
obj_names <- load(here::here("data", "ELSOC_Long.RData"))

raw <- if ("elsoc_long_2016_2022" %in% obj_names) {
  get("elsoc_long_2016_2022")
} else if (length(obj_names) == 1) {
  get(obj_names[1])
} else {
  stop("Cannot identify ELSOC object. Available: ",
       paste(obj_names, collapse = ", "))
}
dt_raw <- dplyr::as_tibble(raw)
message("Loaded: ", nrow(dt_raw), " obs x ", ncol(dt_raw), " variables")

# ── 2. Sample filters (muestra == 1, tipo_atricion == 1) ───────────────────
dt <- dt_raw

if (FILTER_MUESTRA && "muestra" %in% names(dt)) {
  n_antes <- nrow(dt)
  dt <- dt %>% dplyr::filter(muestra == 1)
  message("Filter muestra==1: ", n_antes, " -> ", nrow(dt), " obs")
} else if (FILTER_MUESTRA) {
  warning("Variable 'muestra' not found. Filter not applied.")
}

if (FILTER_TIPO_ATRICION && "tipo_atricion" %in% names(dt)) {
  n_antes <- nrow(dt)
  dt <- dt %>% dplyr::filter(tipo_atricion == 1)
  message("Filter tipo_atricion==1: ", n_antes, " -> ", nrow(dt), " obs")
} else if (FILTER_TIPO_ATRICION) {
  warning("Variable 'tipo_atricion' not found. Filter not applied.")
}

# Filter to analysis waves
dt <- dt %>%
  dplyr::filter(ola %in% WAVES_RAW) %>%
  dplyr::mutate(
    id  = idencuesta,
    ola_rec = dplyr::case_when(
      ola == WAVES_RAW[1] ~ 1L,
      ola == WAVES_RAW[2] ~ 2L,
      ola == WAVES_RAW[3] ~ 3L
    ),
    wave_year = dplyr::case_when(
      ola_rec == 1 ~ WAVES_YEARS[1],
      ola_rec == 2 ~ WAVES_YEARS[2],
      ola_rec == 3 ~ WAVES_YEARS[3]
    )
  ) %>%
  dplyr::rename(ola_orig = ola, ola = ola_rec)

message("Observations in analysis waves (1,3,6): ", nrow(dt))

# ── 3. Membership indicators ───────────────────────────────────────────
missing_member <- setdiff(MEMBER_ITEMS, names(dt))
if (length(missing_member) > 0)
  stop("Membership items not found: ", paste(missing_member, collapse=", "))

dt <- dt %>%
  dplyr::mutate(
    nhg          = member_binary(c12_01),
    religious    = member_binary(c12_02),
    sport        = member_binary(c12_03),
    charity      = member_binary(c12_04),
    political    = member_binary(c12_05),
    professional = member_binary(c12_06),
    union        = member_binary(c12_07),
    student      = member_binary(c12_08)
  )

# Check whether c12_09 (9th domain) exists in the dataset
if (MEMBER_ITEMS_9TH %in% names(dt)) {
  dt <- dt %>% dplyr::mutate(otra = member_binary(c12_09))
  if (!"otra" %in% DOMAIN_VARS) DOMAIN_VARS <<- c(DOMAIN_VARS, "otra")
  N_DOMAINS <<- 9L
  message("c12_09 found -- added as 9th domain.")
} else {
  message("c12_09 not found. Using 8 domains.")
}

# ── 4. Outcome variables ──────────────────────────────────────────────
for (v in c(TRUST_GEN_VAR, TRUST_NH_VAR)) {
  if (!v %in% names(dt)) stop("Variable not found: ", v,
                               "\nExpected t01 for neighbourhood trust -- check codebook.")
}

dt <- dt %>%
  dplyr::mutate(
    trust    = code_gen_trust(.data[[TRUST_GEN_VAR]]),
    trust_nh = code_nh_trust(.data[[TRUST_NH_VAR]])
  )

message("Trust generalizada: ", round(mean(dt$trust, na.rm=TRUE), 3),
        " | Trust vecinal: ", round(mean(dt$trust_nh, na.rm=TRUE), 3))

# ── 5. Covariates (replicating estimation.do) ─────────────────────────────
dt <- dt %>%
  dplyr::mutate(
    # Sex: m0_sexo (1=male, 2=female in ELSOC)
    woman     = dplyr::case_when(m0_sexo == 1 ~ 0L, m0_sexo == 2 ~ 1L, TRUE ~ NA_integer_),
    edad      = to_na(m0_edad),
    education = code_education_5(.data[[EDU_VAR]]),  # 5 niveles como Stata
    edu_bin   = code_edu_binary(.data[[EDU_VAR]]),

    # employed = 1 if m02 <= 2
    employed  = if (EMPLOY_VAR %in% names(.)) {
      code_employment(.data[[EMPLOY_VAR]])
    } else {
      warning("'", EMPLOY_VAR, "' not found."); NA_integer_
    },

    # swb = 1 if s01 > 4
    swb = if (SWB_VAR %in% names(.)) {
      code_swb(.data[[SWB_VAR]])
    } else {
      warning("'", SWB_VAR, "' not found."); NA_integer_
    },

    # couple = 1 if m36 < 4
    couple = if (COUPLE_VAR %in% names(.)) {
      code_couple(.data[[COUPLE_VAR]])
    } else {
      warning("'", COUPLE_VAR, "' not found."); NA_integer_
    },

    # tinst: institutional trust (raw; control in Stata models)
    tinst = if (TINST_VAR %in% names(.)) {
      code_tinst(.data[[TINST_VAR]])
    } else {
      warning("'", TINST_VAR, "' not found."); NA_real_
    }
  )

# ── 6. Derived membership variables ───────────────────────────────────────
dt <- dt %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    membership_count = sum(dplyr::c_across(dplyr::all_of(DOMAIN_VARS)), na.rm = TRUE),
    domain_diversity = membership_count / N_DOMAINS,
    multi_domain     = as.integer(membership_count >= 2L)
  ) %>%
  dplyr::ungroup()

# ── 7. Panel balancing ───────────────────────────────────────────────────────────────
# (a) complete.cases() on model variables (equiv. probit e(sample))
# (b) keep only ids with exactly 3 waves

# Step (a): variables used by Stata to define the analytic sample
stata_model_vars <- c("trust","swb","employed","couple","tinst","woman","edad")
stata_model_vars <- intersect(stata_model_vars, names(dt))
stata_model_vars <- stata_model_vars[sapply(stata_model_vars,
                    function(v) any(!is.na(dt[[v]])))]

# Observations with complete data on all model variables
complete_mask <- complete.cases(dt[, stata_model_vars])
dt$in_model   <- complete_mask

message("Step (a) -- complete cases on model vars: ",
        sum(complete_mask), " de ", nrow(dt), " obs (",
        round(mean(complete_mask)*100, 1), "%)")

# Step (b): keep only ids observed in all 3 waves
panel_counts <- dt %>%
  dplyr::filter(in_model) %>%
  dplyr::count(id)

balanced_ids <- panel_counts %>% dplyr::filter(n == 3) %>% dplyr::pull(id)

dt$in_balanced <- dt$id %in% balanced_ids & dt$in_model

n_balanced <- length(balanced_ids)
n_total    <- length(unique(dt$id))
message("Step (b) -- Balanced panel: n = ", n_balanced, " individuals (",
        round(n_balanced / n_total * 100, 1), "% retention)")

# ── 8. Descriptive statistics ──────────────────────────────────────────────
summary_fun <- function(d) {
  d %>%
    dplyr::group_by(ola) %>%
    dplyr::summarise(
      N                  = dplyr::n(),
      age_mean           = mean(edad,             na.rm = TRUE),
      age_sd             = sd(edad,              na.rm = TRUE),
      pct_woman          = mean(woman,            na.rm = TRUE),
      pct_edu_secondary  = mean(edu_bin,          na.rm = TRUE),
      pct_employed       = mean(employed,         na.rm = TRUE),
      pct_swb            = mean(swb,              na.rm = TRUE),
      pct_couple         = mean(couple,           na.rm = TRUE),
      generalized_trust  = mean(trust,            na.rm = TRUE),
      neighborhood_trust = mean(trust_nh,         na.rm = TRUE),
      membership_mean    = mean(membership_count, na.rm = TRUE),
      membership_sd      = sd(membership_count,  na.rm = TRUE),
      pct_multi_domain   = mean(multi_domain,     na.rm = TRUE),
      domain_diversity   = mean(domain_diversity, na.rm = TRUE),
      .groups = "drop"
    )
}

stats_balanced <- dt %>% dplyr::filter(in_balanced) %>% summary_fun()
stats_full     <- dt %>% summary_fun()

readr::write_csv(stats_balanced, here::here("output", "summary_stats_balanced.csv"))
readr::write_csv(stats_full,     here::here("output", "summary_stats_full.csv"))

cat("\n--- Summary statistics (balanced panel) ---\n")
print(stats_balanced)

# ── 9. Membership by domain and wave ───────────────────────────────────────
domain_tbl <- dt %>%
  dplyr::filter(in_balanced) %>%
  dplyr::group_by(ola) %>%
  dplyr::summarise(
    N = dplyr::n(),
    dplyr::across(dplyr::all_of(DOMAIN_VARS), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(dplyr::all_of(DOMAIN_VARS),
                      names_to = "domain", values_to = "pct_member") %>%
  dplyr::arrange(ola, dplyr::desc(pct_member))

readr::write_csv(domain_tbl, here::here("output", "membership_by_domain.csv"))

# ── 10. Distribution checks ─────────────────────────────────────────────────
cat("\n--- Check: Generalized trust (c02), balanced panel, wave 1 ---\n")
dt %>% dplyr::filter(in_balanced, ola == 1) %>%
  dplyr::count(.data[[TRUST_GEN_VAR]], trust) %>% print()

cat("\n--- Check: Neighbourhood trust (", TRUST_NH_VAR, "), balanced panel, wave 1 ---\n")
dt %>% dplyr::filter(in_balanced, ola == 1) %>%
  dplyr::count(.data[[TRUST_NH_VAR]], trust_nh) %>% print()

cat("\n--- Check: swb, balanced panel, wave 1 ---\n")
dt %>% dplyr::filter(in_balanced, ola == 1) %>%
  dplyr::count(.data[[SWB_VAR]], swb) %>% print()

cat("\n--- Check: Membership (", MEMBER_CODE_LOGIC, "), wave 1 ---\n")
cat("Share with >= 1 membership:",
    round(mean(dt$membership_count[dt$in_balanced & dt$ola==1] > 0, na.rm=TRUE), 3), "\n")

# ── 11. Save analysis dataset ──────────────────────────────────────────────
dt_analysis <- dt %>%
  dplyr::filter(in_balanced) %>%
  dplyr::select(
    id, ola, wave_year, trust, trust_nh, woman, edad, education, edu_bin,
    employed, swb, couple, tinst,
    dplyr::all_of(DOMAIN_VARS),
    membership_count, domain_diversity, multi_domain
  )

saveRDS(dt_analysis, here::here("data", "dt_analysis.rds"))
message("\n[01_descriptive_stats.R] done.")
message("  Balanced panel: n = ", length(balanced_ids),
        " individuals | ", nrow(dt_analysis), " person-wave obs")
message("  Saved: data/dt_analysis.rds")
message("  If trust_nh is 0% or 100%, check that t01 exists in the dataset.")
message("  If employed is all NA, check that m02 exists.")
message("  MEMBER_CODE_LOGIC = '", MEMBER_CODE_LOGIC, "'")
