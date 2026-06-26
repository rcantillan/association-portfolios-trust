# ── 03_trust_models.R — RE probit trust models ───────────────────────────────
#
# Canonical STATE_MAP (from 00_setup.R):
#   α (alpha) = Isolate/Apathetic  | clase=3 | ref=ib3
#   β (beta)  = Closed/Clustering  | clase=1 | ref=ib1
#   γ (gamma) = Bridging/Broker    | clase=2 | ref=ib2 [REFERENCE in estimation.do]
#
# Methodological notes:
#   nAGQ=12: matches Stata xtprobit GHQ-12 quadrature (Laplace default is biased).
#   Cluster-robust SEs via sandwich::vcovCL (requires merDeriv for estfun.glmerMod).
#   marginaleffects_safe=FALSE: silences warning about fixed-effect VCOV in GLMM,
#     consistent with Stata margins post-xtprobit.
#
# Outputs:
#   output/pre_model_diagnostics.txt
#   output/table4_reprobit_ref_bridging.csv
#   output/table5_reprobit_ref_closed.csv
#   output/tableA1_probit_pool_ref_closed.csv
#   output/tableA2_reverse_causality.csv
#   output/fig3_predicted_probs.png / .pdf
#   output/sensitivity_confgral.csv
#   output/ame_all_models.csv
#   output/trust_models_summaries.txt
#   data/trust_model_objects.rds

options(marginaleffects_safe = FALSE)

source(here::here("code", "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(lme4); library(marginaleffects)
  library(sandwich); library(ggplot2); library(broom.mixed)
  library(stringr)
})

# merDeriv provides estfun() and bread() for merMod objects, required for
# sandwich::vcovCL to work with glmer.
if (!requireNamespace("merDeriv", quietly = TRUE))
  install.packages("merDeriv", quiet = TRUE)
suppressPackageStartupMessages(library(merDeriv))

dir.create(here::here("output"), showWarnings = FALSE)
dir.create(here::here("data"),   showWarnings = FALSE)

# nAGQ=12 matches Stata xtprobit GHQ-12 (nAGQ=9 is faster but slightly less exact)
NAGQ <- 12L

# ── 1. Load and merge ─────────────────────────────────────────────────────────
stop_if_missing(c(
  here::here("data", "dt_states_cov.rds"),
  here::here("data", "dt_analysis.rds")
))

dt_states <- readRDS(here::here("data", "dt_states_cov.rds")) %>% as_tibble()
dt_anal   <- readRDS(here::here("data", "dt_analysis.rds"))   %>% as_tibble()

message("dt_states_cov: ", nrow(dt_states), " obs | ",
        length(unique(dt_states$idencuesta)), " ids")
message("dt_analysis:   ", nrow(dt_anal), " obs | ",
        length(unique(dt_anal$id)), " ids")

dt_anal_join <- dt_anal %>%
  dplyr::rename(idencuesta = id) %>%
  dplyr::select(
    idencuesta, ola,
    trust, trust_nh,
    edad, woman, education, edu_bin,
    employed, swb, couple, tinst,
    membership_count, domain_diversity, multi_domain
  )

dt_states    <- dt_states    %>% dplyr::mutate(idencuesta = as.character(idencuesta))
dt_anal_join <- dt_anal_join %>% dplyr::mutate(idencuesta = as.character(idencuesta))

# Only bring columns from dt_anal not already in dt_states (avoid .x/.y conflicts)
cols_new <- setdiff(names(dt_anal_join),
                    c("idencuesta", "ola", setdiff(names(dt_states), c("idencuesta","ola"))))
dt <- dplyr::left_join(
  dt_states,
  dt_anal_join %>% dplyr::select(dplyr::all_of(c("idencuesta", "ola", cols_new))),
  by = c("idencuesta", "ola")
)

n_matched   <- sum(!is.na(dt$trust))
pct_matched <- round(n_matched / nrow(dt) * 100, 1)
message("Match: ", n_matched, "/", nrow(dt), " obs (", pct_matched, "%)")
if (pct_matched < 80)
  warning("< 80% obs have trust. Check filters in 01 and 02.")

# ── 2. Derived variables — canonical STATE_MAP ────────────────────────────────
dt <- dt %>%
  dplyr::mutate(
    clase = dplyr::case_when(
      position == "alpha" ~ 3L,   # α Isolate   clase=3 ib3
      position == "beta"  ~ 1L,   # β Closed    clase=1 ib1
      position == "gamma" ~ 2L,   # γ Bridging  clase=2 ib2 [reference]
      TRUE                ~ NA_integer_
    ),
    clase_strict = dplyr::case_when(
      position_strict == "alpha" ~ 3L,
      position_strict == "beta"  ~ 1L,
      position_strict == "gamma" ~ 2L,
      TRUE                       ~ NA_integer_
    ),
    clase_label = dplyr::recode(as.character(clase), !!!STATE_LABELS_FIG)
  )

# CONTROLS_TRUST replicates estimation.do: swb, employed, couple.
# Education excluded to avoid collinearity with portfolio class.
CONTROLS_TRUST <- c("swb", "employed", "couple")

ctrl_use <- CONTROLS_TRUST[
  CONTROLS_TRUST %in% names(dt) &
  sapply(CONTROLS_TRUST, function(v) sum(!is.na(dt[[v]])) > 0)
]
message("Controls available: ", paste(ctrl_use, collapse = ", "))

# ── 3. Pre-model diagnostics ──────────────────────────────────────────────────
# Review output/pre_model_diagnostics.txt before interpreting results.
# Stata targets: N = 3,891 | Pattern: γ(Bridging) > α(Isolate) > β(Closed)
cat("\n=== PRE-MODEL DIAGNOSTICS ===\n")

diag_lines <- c(
  "=== PRE-MODEL DIAGNOSTICS (03_trust_models.R) ===",
  paste0("Fecha: ", Sys.time()),
  paste0("nAGQ = ", NAGQ, " | Stata GHQ = 12"),
  ""
)

# 3a. Sample size and loss breakdown
stata_model_vars_check <- intersect(c("trust","clase",ctrl_use,"ola"), names(dt))
d_check <- dt %>%
  dplyr::filter(!is.na(clase)) %>%
  tidyr::drop_na(dplyr::all_of(stata_model_vars_check))

cat("N available for main model (trust + controls): ", nrow(d_check), "\n")
cat("N target (Stata):                                 3,891\n")
cat("Difference:                                      ", nrow(d_check) - 3891, "\n\n")

# NA breakdown by variable
cat("=== NA counts by variable (in dt with non-NA clase) ===\n")
dt_con_clase <- dt %>% dplyr::filter(!is.na(clase))
cat("  N total with non-NA clase: ", nrow(dt_con_clase), "\n")
for (v in stata_model_vars_check) {
  n_miss <- sum(is.na(dt_con_clase[[v]]))
  cat("  NA en", sprintf("%-15s", v), ":", n_miss, "\n")
}
cat("\n")
cat("  If difference vs Stata (3,891) > 200:\n")
cat("  -> Check muestra==1 and tipo_atricion==1 filters in 01_descriptive_stats.R\n")
cat("  -> Check that dt_states + dt_analysis merge does not drop obs\n\n")

diag_lines <- c(diag_lines,
  paste0("N disponible R (trust+controles): ", nrow(d_check)),
  paste0("N total con clase no-NA: ", nrow(dt_con_clase)),
  paste0("N target Stata: 3,891"),
  paste0("Difference: ", nrow(d_check) - 3891),
  "  If |diff| > 200, check muestra/tipo_atricion filters in 01.",
  ""
)

# 3b. Mean trust by class — expected pattern: γ > α > β
cat("=== Trust by class (expected: Broker[γ] > Apathetic[α] > Closed[β]) ===\n")
trust_by_clase <- dt %>%
  dplyr::filter(!is.na(trust), !is.na(clase)) %>%
  dplyr::group_by(position, clase, clase_label) %>%
  dplyr::summarise(
    n        = dplyr::n(),
    trust    = round(mean(trust,    na.rm = TRUE), 3),
    trust_nh = round(mean(trust_nh, na.rm = TRUE), 3),
    .groups  = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(trust))
print(trust_by_clase)

diag_lines <- c(diag_lines,
  "--- Trust por clase (esperado: γ > α > β) ---",
  capture.output(print(trust_by_clase)), ""
)

# 3c. Class distribution by wave
cat("\n=== Class distribution by wave ===\n")
clase_dist <- dt %>%
  dplyr::filter(!is.na(clase)) %>%
  dplyr::count(ola, clase_label) %>%
  tidyr::pivot_wider(names_from = ola, values_from = n, names_prefix = "ola_")
print(clase_dist)

# 3d. Trust variable coding check
cat("\n=== trust (expected: binary 0/1) ===\n")
print(table(dt$trust,    useNA = "ifany"))
print(table(dt$trust_nh, useNA = "ifany"))

diag_lines <- c(diag_lines,
  "--- Trust distribution ---",
  capture.output(print(table(dt$trust,    useNA="ifany"))),
  capture.output(print(table(dt$trust_nh, useNA="ifany"))), ""
)

writeLines(diag_lines, here::here("output", "pre_model_diagnostics.txt"))
message("  Saved: output/pre_model_diagnostics.txt")

# ── 4. Helpers ─────────────────────────────────────────────────────────────────

# Returns cluster formula (~idencuesta_f) for avg_slopes(vcov=),
# or NULL if merDeriv is unavailable (fallback to model vcov).
# Uses sandwich::vcovCL with HC1 correction (clubSandwich fails on S4 glmerMod).
get_vcov_formula <- function(model, cluster_var = "idencuesta_f") {
  # Verify merDeriv provides estfun for this model
  has_estfun <- tryCatch({
    sf <- sandwich::estfun(model)
    is.matrix(sf) && nrow(sf) > 0
  }, error = function(e) FALSE)

  if (!has_estfun) {
    message("  estfun() failed for glmerMod -- merDeriv not available or incompatible")
    message("  -> fallback to model vcov (non-clustered SEs)")
    return(NULL)
  }

  # Test that vcovCL works before returning the formula
  vcv_test <- tryCatch(
    sandwich::vcovCL(model,
                     cluster = as.formula(paste0("~", cluster_var)),
                     type    = "HC1"),
    error   = function(e) { message("  vcovCL error: ", e$message); NULL },
    warning = function(w) { message("  vcovCL warning: ", w$message); NULL }
  )

  if (is.null(vcv_test)) {
    message("  vcovCL failed -> fallback to model vcov")
    return(NULL)
  }

  # Validar PD
  ev <- tryCatch(
    eigen(vcv_test, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(ev) || any(ev < -1e-8)) {
    message("  vcov HC1 no-PD -> fallback")
    return(NULL)
  }

  message("  vcovCL HC1 OK -- using cluster-robust SEs")
  as.formula(paste0("~", cluster_var))
}

# RE probit with nAGQ (equivalent to Stata xtprobit GHQ)
fit_reprobit <- function(outcome, data, ref_clase, with_controls, label) {

  ref_label <- dplyr::filter(STATE_MAP, clase == ref_clase) %>%
    dplyr::pull(label_long)
  message("\n  Fitting: ", label,
          " | outcome=",  outcome,
          " | ref=clase", ref_clase, " (", ref_label, ")",
          " | nAGQ=", NAGQ)

  d <- data %>%
    dplyr::filter(!is.na(clase), !is.na(.data[[outcome]])) %>%
    dplyr::mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola),
      clase_f      = relevel(as.factor(clase), ref = as.character(ref_clase))
    )

  ctrl <- if (with_controls) {
    ctrl_use[ctrl_use %in% names(d) &
             sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  } else character(0)

  d <- tidyr::drop_na(
    d, dplyr::all_of(intersect(c(outcome, "clase_f", "ola_f", ctrl), names(d)))
  )

  if (nrow(d) < 50) { warning("n muy bajo en ", label, ": ", nrow(d)); return(NULL) }
  message("    N obs=", nrow(d), " | N ids=", length(unique(d$idencuesta_f)))

  rhs <- paste(c("clase_f", ctrl, "ola_f", "(1 | idencuesta_f)"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  # nAGQ >= 2 with 1 RE is equivalent to Stata GHQ.
  # Falls back to Laplace (nAGQ=1) on convergence error.
  m <- tryCatch(
    lme4::glmer(
      fml, data = d,
      family  = binomial(link = "probit"),
      nAGQ    = NAGQ,
      control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl   = list(maxfun = 5e5)
      )
    ),
    error = function(e) {
      message("  glmer nAGQ=", NAGQ, " failed. Retrying with Laplace (nAGQ=1).")
      lme4::glmer(
        fml, data = d,
        family  = binomial(link = "probit"),
        nAGQ    = 1L,
        control = lme4::glmerControl(optimizer = "bobyqa",
                                     optCtrl   = list(maxfun = 5e5))
      )
    }
  )

  sigma_u <- round(sqrt(as.numeric(VarCorr(m)$idencuesta_f)), 4)
  message("    sigma_u=", sigma_u, " (Stata typical: 1.5-2.5 for RE probit on trust)")

  vcov_fml  <- get_vcov_formula(m)
  vcov_type <- ifelse(is.null(vcov_fml), "model_vcov (fallback)", "HC1_cluster")

  # vcov = formula ~cluster -> sandwich::vcovCL (cluster-robust HC1)
  # vcov = NULL             -> marginaleffects uses vcov(model) (fallback)
  ame <- marginaleffects::avg_slopes(
    m, variables = "clase_f", type = "response",
    vcov = if (is.null(vcov_fml)) TRUE else vcov_fml
  ) %>%
    as_tibble() %>%
    dplyr::mutate(
      outcome       = outcome,
      ref_clase     = ref_clase,
      ref_label     = ref_label,
      with_controls = with_controls,
      spec          = label,
      n_obs         = nrow(d),
      n_id          = length(unique(d$idencuesta_f)),
      vcov_type     = vcov_type,
      sigma_u       = sigma_u
    )

  # Label contrasts with canonical Greek symbols
  ame <- ame %>%
    dplyr::mutate(
      contrast_label = stringr::str_replace_all(
        contrast,
        c(
          "mean\\(1\\)" = "\u03b2 Closed",
          "mean\\(2\\)" = "\u03b3 Bridging",
          "mean\\(3\\)" = "\u03b1 Isolate",
          "^1 - 2$"     = "\u03b2 vs \u03b3",
          "^3 - 2$"     = "\u03b1 vs \u03b3",
          "^2 - 1$"     = "\u03b3 vs \u03b2",
          "^3 - 1$"     = "\u03b1 vs \u03b2",
          "^1 - 3$"     = "\u03b2 vs \u03b1",
          "^2 - 3$"     = "\u03b3 vs \u03b1"
        )
      )
    )

  # Predicted probabilities by class x wave for Figure 3
  pred <- tryCatch(
    marginaleffects::predictions(
      m,
      newdata = marginaleffects::datagrid(
        clase_f   = levels(d$clase_f),
        ola_f     = levels(d$ola_f),
        grid_type = "mean_or_mode"
      ),
      type = "response"
    ) %>%
      as_tibble() %>%
      dplyr::mutate(
        outcome     = outcome,
        spec        = label,
        clase_int   = as.integer(as.character(clase_f)),
        clase_label = dplyr::recode(as.character(clase_int), !!!STATE_LABELS_FIG)
      ),
    error = function(e) {
      warning("predictions() failed in ", label, ": ", e$message)
      NULL
    }
  )

  list(model = m, ame = ame, pred = pred,
       formula   = deparse(fml), n = nrow(d),
       label     = label, ref_label = ref_label,
       vcov_type = vcov_type, sigma_u = sigma_u)
}

# Pooled probit helper (Table A1)
fit_probit_pool <- function(outcome, data, ref_clase, with_controls, label) {

  ref_label <- dplyr::filter(STATE_MAP, clase == ref_clase) %>%
    dplyr::pull(label_long)

  d <- data %>%
    dplyr::filter(!is.na(clase), !is.na(.data[[outcome]])) %>%
    dplyr::mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola),
      clase_f      = relevel(as.factor(clase), ref = as.character(ref_clase))
    )

  ctrl <- if (with_controls) {
    ctrl_use[ctrl_use %in% names(d) &
             sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  } else character(0)

  d <- tidyr::drop_na(
    d, dplyr::all_of(intersect(c(outcome, "clase_f", "ola_f", ctrl), names(d)))
  )

  rhs <- paste(c("clase_f", ctrl, "ola_f"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  m   <- glm(fml, data = d, family = binomial(link = "probit"))
  vcv <- tryCatch(sandwich::vcovCL(m, cluster = ~idencuesta_f),
                  error = function(e) NULL)

  ame <- marginaleffects::avg_slopes(
    m, variables = "clase_f", type = "response", vcov = vcv
  ) %>%
    as_tibble() %>%
    dplyr::mutate(
      outcome = outcome, ref_clase = ref_clase, ref_label = ref_label,
      with_controls = with_controls, spec = label, n_obs = nrow(d)
    )

  list(model = m, ame = ame, formula = deparse(fml), label = label)
}

# Reverse causality helper: trust -> clase
fit_reverse <- function(trust_var, data, with_controls, label) {

  d <- data %>%
    dplyr::filter(!is.na(clase), !is.na(.data[[trust_var]])) %>%
    dplyr::mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola)
    )

  ctrl <- if (with_controls) {
    ctrl_use[ctrl_use %in% names(d) &
             sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  } else character(0)

  d <- tidyr::drop_na(
    d, dplyr::all_of(intersect(c("clase", trust_var, ctrl), names(d)))
  )

  rhs <- paste(c(trust_var, ctrl, "ola_f", "(1 | idencuesta_f)"), collapse = " + ")
  fml <- as.formula(paste("clase ~", rhs))

  m <- lme4::lmer(fml, data = d,
                  control = lme4::lmerControl(optimizer = "bobyqa"))

  vcv <- tryCatch({
    v  <- sandwich::vcovCL(m, cluster = d$idencuesta_f, type = "HC1")
    ev <- eigen(v, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < -1e-8)) NULL else v
  }, error = function(e) {
    message("  vcovCL error en fit_reverse: ", e$message, " -> NULL")
    NULL
  })

  slopes <- tryCatch(
    marginaleffects::avg_slopes(m, variables = trust_var, vcov = vcv) %>%
      as_tibble() %>%
      dplyr::mutate(
        predictor     = trust_var,
        with_controls = with_controls,
        spec          = label,
        n_obs         = nrow(d)
      ),
    error = function(e) NULL
  )

  list(model = m, slopes = slopes, formula = deparse(fml), label = label)
}

# ── 5. Table 4 — RE probit ref=γ Bridging (ib2.clase) ────────────────────────
message("\n=== TABLE 4: RE probit | ref=\u03b3 Bridging | nAGQ=", NAGQ, " ===")
# Stata targets: β=-0.071***, α=-0.040* (M1); β=-0.060***, α=-0.037* (M2);
#                β=-0.046(ns), α=-0.091*** (M3); β=-0.073***, α=-0.092*** (M4)

t4_m1 <- fit_reprobit("trust",    dt, ref_clase=2, with_controls=FALSE, label="T4_M1")
t4_m2 <- fit_reprobit("trust",    dt, ref_clase=2, with_controls=TRUE,  label="T4_M2")
t4_m3 <- fit_reprobit("trust_nh", dt, ref_clase=2, with_controls=FALSE, label="T4_M3")
t4_m4 <- fit_reprobit("trust_nh", dt, ref_clase=2, with_controls=TRUE,  label="T4_M4")

table4 <- dplyr::bind_rows(t4_m1$ame, t4_m2$ame, t4_m3$ame, t4_m4$ame)
readr::write_csv(table4, here::here("output", "table4_reprobit_ref_bridging.csv"))

cat("\n--- Table 4 (ref=\u03b3 Bridging; negativo = menos trust que Bridging) ---\n")
table4 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value < .001 ~ "***", p.value < .01 ~ "**",
    p.value < .05  ~ "*",   p.value < .10  ~ ".", TRUE ~ ""
  )) %>%
  dplyr::select(spec, outcome, contrast_label, estimate, std.error, p.value, sig,
                n_obs, vcov_type, sigma_u) %>%
  print()

# Direct comparison with Stata targets
cat("\n--- Comparison vs Stata (N=3,891, Table 4) ---\n")
# Stata targets for Table 4 (ref=γ Bridging)
stata_ref <- tibble::tibble(
  spec           = c("T4_M1","T4_M1","T4_M2","T4_M2","T4_M3","T4_M3","T4_M4","T4_M4"),
  contrast_label = rep(c("\u03b2 vs \u03b3", "\u03b1 vs \u03b3"), 4),
  stata_est      = c(-0.071,-0.040, -0.060,-0.037, -0.046,-0.091, -0.073,-0.092),
  stata_se       = c(0.021, 0.020,  0.020, 0.020,  0.032, 0.030,  0.032, 0.031),
  stata_sig      = c("***","*",    "***","*",     "",    "***",   "***","***")
)
table4 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::select(spec, contrast_label, estimate, std.error, p.value) %>%
  dplyr::mutate(
    sig   = dplyr::case_when(p.value<.001~"***",p.value<.01~"**",
                             p.value<.05~"*",p.value<.10~".",TRUE~""),
    R_est = round(estimate,  3),
    R_se  = round(std.error, 3)
  ) %>%
  dplyr::left_join(stata_ref, by = c("spec", "contrast_label")) %>%
  dplyr::select(spec, contrast_label, R_est, R_se, sig, stata_est, stata_se, stata_sig) %>%
  print()

# ── 6. Table 5 — RE probit ref=β Closed (ib1.clase) ──────────────────────────
message("\n=== TABLE 5: RE probit | ref=\u03b2 Closed | nAGQ=", NAGQ, " ===")

t5_m5 <- fit_reprobit("trust",    dt, ref_clase=1, with_controls=FALSE, label="T5_M5")
t5_m6 <- fit_reprobit("trust",    dt, ref_clase=1, with_controls=TRUE,  label="T5_M6")
t5_m7 <- fit_reprobit("trust_nh", dt, ref_clase=1, with_controls=FALSE, label="T5_M7")
t5_m8 <- fit_reprobit("trust_nh", dt, ref_clase=1, with_controls=TRUE,  label="T5_M8")

table5 <- dplyr::bind_rows(t5_m5$ame, t5_m6$ame, t5_m7$ame, t5_m8$ame)
readr::write_csv(table5, here::here("output", "table5_reprobit_ref_closed.csv"))

cat("\n--- Table 5 (ref=\u03b2 Closed) ---\n")
table5 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",p.value<.10~".",TRUE~""
  )) %>%
  dplyr::select(spec, outcome, contrast_label, estimate, std.error, p.value, sig) %>%
  print()

# ── 7. Table A1 — Pooled probit ref=β Closed ─────────────────────────────────
message("\n=== TABLE A1: Pooled probit | ref=\u03b2 Closed ===")

a1_m1 <- fit_probit_pool("trust",    dt, ref_clase=1, with_controls=FALSE, label="A1_M1")
a1_m2 <- fit_probit_pool("trust",    dt, ref_clase=1, with_controls=TRUE,  label="A1_M2")
a1_m3 <- fit_probit_pool("trust_nh", dt, ref_clase=1, with_controls=FALSE, label="A1_M3")
a1_m4 <- fit_probit_pool("trust_nh", dt, ref_clase=1, with_controls=TRUE,  label="A1_M4")

tableA1 <- dplyr::bind_rows(a1_m1$ame, a1_m2$ame, a1_m3$ame, a1_m4$ame)
readr::write_csv(tableA1, here::here("output", "tableA1_probit_pool_ref_closed.csv"))

cat("\n--- Table A1 ---\n")
tableA1 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",TRUE~""
  )) %>%
  dplyr::select(spec, outcome, contrast, estimate, std.error, p.value, sig) %>%
  print()

# ── 8. Sensitivity: conf_gral alternative outcome ─────────────────────────────
message("\n=== SENSIBILIDAD: conf_gral ===")
if ("conf_gral" %in% names(dt) && sum(!is.na(dt$conf_gral)) > 100) {
  s_b <- fit_reprobit("conf_gral", dt, ref_clase=2, with_controls=TRUE, label="Sens_Bridging")
  s_c <- fit_reprobit("conf_gral", dt, ref_clase=1, with_controls=TRUE, label="Sens_Closed")
  if (!is.null(s_b) && !is.null(s_c)) {
    dplyr::bind_rows(s_b$ame, s_c$ame) %>%
      readr::write_csv(here::here("output","sensitivity_confgral.csv"))
  }
} else message("  conf_gral no disponible.")

# ── 9. Table A2 — Reverse causality (trust -> clase) ─────────────────────────
message("\n=== TABLE A2: Causalidad reversa (trust -> clase) ===")

ra1 <- fit_reverse("trust",    dt, with_controls=FALSE, label="A2_M1")
ra2 <- fit_reverse("trust",    dt, with_controls=TRUE,  label="A2_M2")
ra3 <- fit_reverse("trust_nh", dt, with_controls=FALSE, label="A2_M3")
ra4 <- fit_reverse("trust_nh", dt, with_controls=TRUE,  label="A2_M4")

tableA2 <- dplyr::bind_rows(ra1$slopes, ra2$slopes, ra3$slopes, ra4$slopes)
readr::write_csv(tableA2, here::here("output","tableA2_reverse_causality.csv"))

cat("\n--- Table A2 ---\n")
tableA2 %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",p.value<.10~".",TRUE~""
  )) %>%
  dplyr::select(spec, predictor, estimate, std.error, p.value, sig) %>%
  print()

# ── 10. Figure 3 — Pr(trust=1) by class and wave ─────────────────────────────
message("\n=== FIGURE 3 ===")

prep_pred <- function(pred_tbl, outcome_label) {
  if (is.null(pred_tbl)) return(NULL)
  pred_tbl %>%
    dplyr::mutate(
      wave_year     = dplyr::case_when(
        ola_f == "1" ~ WAVES_YEARS[1],
        ola_f == "2" ~ WAVES_YEARS[2],
        ola_f == "3" ~ WAVES_YEARS[3],
        TRUE         ~ as.integer(as.character(ola_f))
      ),
      outcome_label = outcome_label
    )
}

preds <- dplyr::bind_rows(
  prep_pred(t4_m2$pred, "Generalized trust (c02)"),
  prep_pred(t4_m4$pred, "Trust in neighbors (t01)")
)

if (!is.null(preds) && nrow(preds) > 0) {

  if (!"clase_label" %in% names(preds)) {
    preds <- preds %>%
      dplyr::mutate(
        clase_int   = as.integer(as.character(clase_f)),
        clase_label = dplyr::recode(as.character(clase_int), !!!STATE_LABELS_FIG)
      )
  }

  # x-axis order: β (Clustering) → γ (Bridging) → α (Isolation)
  preds <- preds %>%
    dplyr::mutate(
      clase_order = dplyr::case_when(
        clase_int == 1L ~ 1L,
        clase_int == 2L ~ 2L,
        clase_int == 3L ~ 3L,
        TRUE ~ NA_integer_
      ),
      wave_label = as.character(wave_year),
      # Short x-axis labels using Greek notation
      clase_label_short = dplyr::case_when(
        clase_int == 1L ~ "β  Clustering",
        clase_int == 2L ~ "γ  Bridging",
        clase_int == 3L ~ "α  Isolation",
        TRUE ~ clase_label
      ),
      # Clean facet labels
      outcome_clean = dplyr::recode(outcome_label,
        "Generalized trust (c02)" = "Generalized trust",
        "Trust in neighbors (t01)" = "Neighborhood trust"
      )
    )

  # Wave colors (NPG, light mint → teal → dark navy = 2016 → 2018 → 2022)
  wave_colors <- c("2016" = "#91D1C2", "2018" = "#4DBBD5", "2022" = "#3C5488")

  fig3 <- ggplot2::ggplot(
    preds,
    ggplot2::aes(
      x    = reorder(clase_label_short, clase_order),
      y    = estimate,
      ymin = conf.low,
      ymax = conf.high,
      fill = wave_label
    )
  ) +
    ggplot2::geom_col(
      position  = ggplot2::position_dodge(width = 0.72),
      width     = 0.62,
      color     = NA
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = conf.low, ymax = conf.high),
      position  = ggplot2::position_dodge(width = 0.72),
      width     = 0.16,
      linewidth = 0.45,
      color     = "grey40"
    ) +
    ggplot2::facet_wrap(~outcome_clean, scales = "fixed") +
    ggplot2::scale_fill_manual(values = wave_colors, name = NULL) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.08))
    ) +
    ggplot2::labs(x = NULL, y = "Predicted probability") +
    theme_ssr() +
    ggplot2::theme(
      legend.position = "top",
      axis.text.x     = ggplot2::element_text(size = 10.5, color = "grey20")
    )

  ggsave_safe("fig3_predicted_probs.png", fig3, width = 8.5, height = 5)
  ggsave_safe("fig3_predicted_probs.pdf", fig3, width = 8.5, height = 5)
  message("  Fig3 saved.")
} else {
  message("  predictions() not available. Fig3 skipped.")
}

# ── 11. Unified AME table ──────────────────────────────────────────────────────
ame_all <- dplyr::bind_rows(
  table4  %>% dplyr::mutate(table = "Table4_ref_Bridging"),
  table5  %>% dplyr::mutate(table = "Table5_ref_Closed"),
  tableA1 %>% dplyr::mutate(table = "TableA1_pool")
) %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",p.value<.10~".",TRUE~""
  )) %>%
  dplyr::arrange(table, outcome, spec, contrast)

readr::write_csv(ame_all, here::here("output","ame_all_models.csv"))

cat("\n=== AMEs significativos (p < .05) ===\n")
ame_all %>%
  dplyr::filter(p.value < .05) %>%
  dplyr::select(table, spec, outcome, contrast_label, estimate, std.error, p.value, sig) %>%
  print(n=40)

# ── 12. Model summaries ────────────────────────────────────────────────────────
mods_to_report <- list(
  "T4_M2 trust ref=\u03b3 +ctrl"    = if (!is.null(t4_m2)) t4_m2$model else NULL,
  "T4_M4 trust_nh ref=\u03b3 +ctrl" = if (!is.null(t4_m4)) t4_m4$model else NULL,
  "T5_M6 trust ref=\u03b2 +ctrl"    = if (!is.null(t5_m6)) t5_m6$model else NULL,
  "T5_M8 trust_nh ref=\u03b2 +ctrl" = if (!is.null(t5_m8)) t5_m8$model else NULL,
  "A2_M2 reverse trust +ctrl"       = if (!is.null(ra2))   ra2$model   else NULL
)

summ_lines <- unlist(lapply(names(mods_to_report), function(nm) {
  m <- mods_to_report[[nm]]
  if (is.null(m)) return(paste0("\n=== ", nm, ": NULL ==="))
  c(paste0("\n\n=== ", nm, " ==="), capture.output(summary(m)))
}))
writeLines(summ_lines, here::here("output","trust_models_summaries.txt"))

# ── 13. Save objects ───────────────────────────────────────────────────────────
saveRDS(
  list(
    dt               = dt,
    ctrl_use         = ctrl_use,
    NAGQ             = NAGQ,
    STATE_MAP        = STATE_MAP,
    STATE_LABELS_FIG = STATE_LABELS_FIG,
    STATE_PALETTE    = STATE_PALETTE,
    t4_m2            = if (!is.null(t4_m2)) t4_m2$model else NULL,
    t4_m4            = if (!is.null(t4_m4)) t4_m4$model else NULL,
    t5_m6            = if (!is.null(t5_m6)) t5_m6$model else NULL,
    t5_m8            = if (!is.null(t5_m8)) t5_m8$model else NULL
  ),
  here::here("data","trust_model_objects.rds")
)

# ── 14. Final summary ──────────────────────────────────────────────────────────
message("\n[03_trust_models.R] done.")
message("  nAGQ=", NAGQ, " | N T4_M2=", if (!is.null(t4_m2)) t4_m2$n else "NULL")
message("  Analysis N target: 1,297 individuals, 3,891 person-waves")
message("  Controls: ", paste(ctrl_use, collapse=", "))
message("")
message("  If results still differ from Stata:")
message("    1. See output/pre_model_diagnostics.txt (N and trust pattern)")
message("    2. Try NAGQ <- 12L for exact GHQ-12 replication")
message("    3. Large N difference -> check filters in 01_descriptive_stats.R")
message("    4. sigma_u << 1.5 -> RE poorly estimated; use nAGQ=12")
message("")
message("  Outputs:")
message("    output/pre_model_diagnostics.txt")
message("    output/table4_reprobit_ref_bridging.csv")
message("    output/table5_reprobit_ref_closed.csv")
message("    output/tableA1_probit_pool_ref_closed.csv")
message("    output/tableA2_reverse_causality.csv")
message("    output/ame_all_models.csv")
message("    output/fig3_predicted_probs.png/.pdf")
