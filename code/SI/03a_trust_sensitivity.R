# ==============================================================================
# 03a_trust_sensitivity.R — Sensitivity Analyses (Supplementary Information)
# ==============================================================================
# Runs five sensitivity checks to accompany the main trust models (03_trust_models.R).
# All outputs go to output/SI_*.
#
# S1 — Correlated Random Effects / Mundlak probit
#        If RE estimates are biased by unobserved heterogeneity, CRE corrects this
#        by adding within-person means of time-varying covariates. If the mean
#        terms are jointly zero, the standard RE is valid.
#
# S2 — Oster (2019) δ (coefficient stability / omitted variable bounds)
#        Asks how strongly a confounder would need to correlate with both clase
#        and trust to fully explain the estimated AME. δ > 1 is the conventional
#        threshold for robustness. Computed via OLS approximation on pooled data.
#
# S3 — Cross-lagged specification (clase_{t-1} → trust_t | trust_{t-1})
#        Stronger test of temporal precedence. Clase at prior wave predicts trust
#        at next wave, controlling for lagged trust. Uses waves 2 and 3 only.
#
# S4 — Entropy-weighted regression (observation weights = p_max)
#        Downweights observations with uncertain latent state assignment.
#        If the effect of clase on trust is driven by confidently classified
#        observations, AMEs should be stable or larger here.
#
# S5 — Soft assignment (p_gamma, p_beta as continuous predictors)
#        Avoids modal classification entirely. The gradient in posterior
#        probabilities directly predicts trust without discretizing positions.
#
# INPUTS:
#   data/dt_states_cov.rds    (from 02_latent_markov.R)
#   data/dt_analysis.rds      (from 01_descriptive_stats.R)
#
# OUTPUTS:
#   output/SI_S1_CRE_mundlak.csv
#   output/SI_S2_oster_delta.csv
#   output/SI_S3_crosslagged.csv
#   output/SI_S4_entropy_weighted.csv
#   output/SI_S5_soft_assignment.csv
#   output/SI_sensitivity_comparison.csv   (all AMEs side by side)
#   output/SI_sensitivity_log.txt
#   output/SI_sensitivity_tables.txt       (LaTeX-ready)
# ==============================================================================

options(marginaleffects_safe = FALSE)

source(here::here("code", "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(lme4); library(marginaleffects); library(clubSandwich)
  library(sandwich); library(stringr); library(broom.mixed)
})

if (!requireNamespace("merDeriv", quietly = TRUE))
  install.packages("merDeriv", quiet = TRUE)
suppressPackageStartupMessages(library(merDeriv))

dir.create(here::here("output"), showWarnings = FALSE)

NAGQ <- 12L   # mismo que 03_trust_models.R

# Log
LOG <- here::here("output", "SI_sensitivity_log.txt")
file.remove(LOG)
lg <- function(...) {
  txt <- paste0(..., collapse = "")
  message(txt)
  write(txt, LOG, append = TRUE)
}

lg("=== 03a_trust_sensitivity.R | ", Sys.time(), " ===")
lg("nAGQ = ", NAGQ)

# ==============================================================================
# 0. CARGA Y PREPARACIÓN (idéntica a 03_trust_models.R)
# ==============================================================================
stop_if_missing(c(
  here::here("data", "dt_states_cov.rds"),
  here::here("data", "dt_analysis.rds")
))

dt_states <- readRDS(here::here("data", "dt_states_cov.rds")) %>% as_tibble()
dt_anal   <- readRDS(here::here("data", "dt_analysis.rds"))   %>% as_tibble()

dt_states    <- dt_states    %>% mutate(idencuesta = as.character(idencuesta))

dt_anal_join <- dt_anal %>%
  rename(idencuesta = id) %>%
  select(idencuesta, ola, trust, trust_nh,
         edad, woman, education, edu_bin,
         employed, swb, couple, tinst,
         membership_count, domain_diversity) %>%
  mutate(idencuesta = as.character(idencuesta))

cols_new <- setdiff(names(dt_anal_join),
                    c("idencuesta", "ola",
                      setdiff(names(dt_states), c("idencuesta","ola"))))
dt <- left_join(
  dt_states,
  dt_anal_join %>% select(all_of(c("idencuesta", "ola", cols_new))),
  by = c("idencuesta", "ola")
)

dt <- dt %>%
  mutate(
    clase = case_when(
      position == "alpha" ~ 3L,
      position == "beta"  ~ 1L,
      position == "gamma" ~ 2L,
      TRUE                ~ NA_integer_
    ),
    clase_label = recode(as.character(clase), !!!STATE_LABELS_FIG)
  )

# Controles base (sin education — ver 03_trust_models.R)
CONTROLS_TRUST <- c("swb", "employed", "couple")

ctrl_use <- CONTROLS_TRUST[
  CONTROLS_TRUST %in% names(dt) &
  sapply(CONTROLS_TRUST, function(v) sum(!is.na(dt[[v]])) > 0)
]
lg("Controles disponibles: ", paste(ctrl_use, collapse = ", "))
lg("N base: ", nrow(dt %>% filter(!is.na(clase), !is.na(trust))), " obs")

# ==============================================================================
# HELPERS
# ==============================================================================

# Etiqueta contrasts con símbolos canónicos
label_contrasts <- function(df) {
  df %>% mutate(
    contrast_label = str_replace_all(
      contrast,
      c("^1 - 2$" = "\u03b2 vs \u03b3", "^3 - 2$" = "\u03b1 vs \u03b3",
        "^2 - 1$" = "\u03b3 vs \u03b2", "^3 - 1$" = "\u03b1 vs \u03b2",
        "^1 - 3$" = "\u03b2 vs \u03b1", "^2 - 3$" = "\u03b3 vs \u03b1",
        "mean\\(1\\)" = "\u03b2", "mean\\(2\\)" = "\u03b3",
        "mean\\(3\\)" = "\u03b1")
    ),
    sig = case_when(
      p.value < .001 ~ "***", p.value < .01 ~ "**",
      p.value < .05  ~ "*",   p.value < .10  ~ ".",
      TRUE           ~ ""
    )
  )
}

# Tabla LaTeX simple
to_latex <- function(df, caption, label) {
  cols <- intersect(c("spec","outcome","contrast_label","estimate",
                      "std.error","p.value","sig"), names(df))
  d <- df %>% select(all_of(cols)) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  header <- paste(cols, collapse = " & ")
  rows <- apply(d, 1, function(r) paste(r, collapse = " & "))
  c(
    paste0("\\begin{table}[H]\\centering"),
    paste0("\\caption{", caption, "}\\label{", label, "}"),
    paste0("\\begin{tabular}{", paste(rep("l", length(cols)), collapse=""), "}"),
    "\\hline\\hline",
    paste0(header, " \\\\"), "\\hline",
    paste0(rows, " \\\\"),
    "\\hline\\hline\\end{tabular}\\end{table}", ""
  )
}

latex_out <- character(0)   # acumulador de tablas LaTeX

# ==============================================================================
# S1 — CORRELATED RANDOM EFFECTS / MUNDLAK PROBIT
# ==============================================================================
lg("\n--- S1: CRE / Mundlak probit ---")

# Variables time-varying para las que agregar la media within-person
tv_vars <- intersect(c("swb", "employed", "couple", "tinst"), names(dt))

dt_cre <- dt %>%
  group_by(idencuesta) %>%
  mutate(across(
    all_of(tv_vars),
    ~ mean(.x, na.rm = TRUE),
    .names = "{.col}_mean"
  )) %>%
  ungroup()

mean_vars   <- paste0(tv_vars, "_mean")
ctrl_cre    <- c(ctrl_use, mean_vars)

fit_cre <- function(outcome, ref_clase, label) {
  # Wooldridge (2010, 2019) CRE / Mundlak pooled probit:
  #   - Pooled probit (glm) + within-person means of time-varying covariates
  #   - Cluster-robust SEs (vcovCL by idencuesta) → valid even if RE structure
  #     is misspecified; gives population-average APEs, not conditional effects
  #   - Wald test of mean terms uses the cluster-robust vcov matrix
  # This is preferred over glmer + RE for APE inference (Wooldridge 2010 §15.8)

  d <- dt_cre %>%
    filter(!is.na(clase), !is.na(.data[[outcome]])) %>%
    mutate(
      ola_f   = as.factor(ola),
      clase_f = relevel(as.factor(clase), ref = as.character(ref_clase))
    )

  ctrl <- ctrl_cre[ctrl_cre %in% names(d) &
                   sapply(ctrl_cre, function(v) sum(!is.na(d[[v]])) > 0)]

  d <- drop_na(d, all_of(intersect(c(outcome, "clase_f", "ola_f", ctrl), names(d))))

  rhs <- paste(c("clase_f", ctrl, "ola_f"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  lg("  Fitting CRE (pooled probit + vcovCL) ", label, " | N=", nrow(d))
  m <- tryCatch(
    glm(fml, data = d, family = binomial("probit")),
    error = function(e) { lg("  CRE glm failed: ", e$message); NULL }
  )
  if (is.null(m)) return(NULL)

  # Cluster-robust vcov (Wooldridge-consistent SEs)
  vcv <- tryCatch(
    sandwich::vcovCL(m, cluster = ~idencuesta),
    error = function(e) { lg("  vcovCL failed: ", e$message); NULL }
  )

  ame <- avg_slopes(m, variables = "clase_f", type = "response",
                    vcov = if (!is.null(vcv)) vcv else TRUE) %>%
    as_tibble() %>%
    mutate(spec = label, outcome = outcome,
           model = "CRE/Mundlak (pooled probit + vcovCL)",
           n_obs = nrow(d))

  # Wald test: mean terms jointly = 0, using cluster-robust vcov
  mean_vars_in <- intersect(mean_vars, names(coef(m)))
  if (length(mean_vars_in) > 0 && !is.null(vcv)) {
    wald_chi2 <- tryCatch({
      vc <- vcv[mean_vars_in, mean_vars_in, drop = FALSE]
      b  <- coef(m)[mean_vars_in]
      as.numeric(t(b) %*% solve(vc) %*% b)
    }, error = function(e) NA_real_)
    wald_p <- 1 - pchisq(wald_chi2, length(mean_vars_in))
    lg("  Wald χ²(", length(mean_vars_in), ") for mean terms = ",
       round(wald_chi2, 2), " (df=", length(mean_vars_in), ")",
       " p=", round(wald_p, 3))
    ame <- ame %>% mutate(
      wald_mean_chi2 = wald_chi2,
      wald_mean_df   = length(mean_vars_in),
      wald_mean_p    = wald_p
    )
  }

  list(model = m, ame = ame)
}

s1_r1 <- fit_cre("trust",    ref_clase = 2L, label = "S1_trust_ref_gamma")
s1_r2 <- fit_cre("trust_nh", ref_clase = 2L, label = "S1_nh_ref_gamma")

s1_ames <- bind_rows(s1_r1$ame, s1_r2$ame) %>%
  filter(!grepl("^ola|_mean", contrast, ignore.case = TRUE)) %>%
  label_contrasts()

write_csv(s1_ames, here::here("output", "SI_S1_CRE_mundlak.csv"))
lg("Saved: SI_S1_CRE_mundlak.csv")

cat("\n--- S1: CRE/Mundlak AMEs (ref = γ) ---\n")
s1_ames %>%
  select(spec, outcome, contrast_label, estimate, std.error, p.value, sig) %>%
  print()

latex_out <- c(latex_out,
  to_latex(s1_ames %>% select(spec, outcome, contrast_label, estimate, std.error, p.value, sig),
           caption = "S1: Correlated Random Effects (Mundlak) pooled probit with cluster-robust standard errors (\\texttt{vcovCL} by respondent). Within-person means of time-varying covariates added as Mundlak corrections. Wald test of joint significance of mean terms reported in text. Reference: $\\gamma$ (Bridging).",
           label   = "tab:SI_S1_CRE"))

# ==============================================================================
# S2 — OSTER (2019) δ
# ==============================================================================
lg("\n--- S2: Oster (2019) delta ---")
#
# Method: OLS pooled regression as approximation (standard in sociology papers).
# R²_max = min(1.3 × R²_full, 1.0)  [Oster's default].
# Computed separately for trust and trust_nh, and for each clase contrast
# (β vs γ, α vs γ) — using indicator variables for each clase.
#
# δ > 1 → confounder would need to explain MORE selection variance than observed
#          controls; conventional threshold for robustness.

oster_rows <- list()

for (outcome_var in c("trust", "trust_nh")) {
  d_ols <- dt %>%
    filter(!is.na(clase), !is.na(.data[[outcome_var]])) %>%
    mutate(
      ola_f   = as.factor(ola),
      beta_d  = as.integer(clase == 1L),   # β indicator
      alpha_d = as.integer(clase == 3L)    # α indicator (γ=2 is ref)
    )

  ctrl_ols <- ctrl_use[ctrl_use %in% names(d_ols) &
                       sapply(ctrl_use, function(v) sum(!is.na(d_ols[[v]])) > 0)]
  d_ols <- drop_na(d_ols, all_of(intersect(
    c(outcome_var, "beta_d", "alpha_d", "ola_f", ctrl_ols), names(d_ols)
  )))

  # Restricted (only clase dummies + wave FE)
  fml_r <- as.formula(paste(outcome_var, "~ beta_d + alpha_d + ola_f"))
  # Full (+ controls)
  fml_f <- as.formula(paste(outcome_var, "~", paste(
    c("beta_d", "alpha_d", ctrl_ols, "ola_f"), collapse = " + "
  )))

  m_r <- lm(fml_r, data = d_ols)
  m_f <- lm(fml_f, data = d_ols)

  r2_r <- summary(m_r)$r.squared
  r2_f <- summary(m_f)$r.squared

  for (contrast_name in c("beta_d", "alpha_d")) {
    beta_r <- coef(m_r)[[contrast_name]]
    beta_f <- coef(m_f)[[contrast_name]]
    delta  <- oster_delta(beta_r, beta_f, r2_r, r2_f)
    label  <- if (contrast_name == "beta_d") "\u03b2 vs \u03b3" else "\u03b1 vs \u03b3"

    oster_rows[[length(oster_rows) + 1]] <- tibble(
      outcome    = outcome_var,
      contrast   = label,
      beta_r     = round(beta_r, 4),
      beta_f     = round(beta_f, 4),
      r2_r       = round(r2_r,   4),
      r2_f       = round(r2_f,   4),
      r2_max     = round(min(1.3 * r2_f, 1.0), 4),
      delta      = round(delta,  3),
      robust     = ifelse(abs(delta) >= 1, "Yes", "No")
    )
    lg("  ", outcome_var, " | ", label,
       " | beta_r=", round(beta_r,3), " beta_f=", round(beta_f,3),
       " | delta=", round(delta,3))
  }
}

s2_tbl <- bind_rows(oster_rows)
write_csv(s2_tbl, here::here("output", "SI_S2_oster_delta.csv"))
lg("Saved: SI_S2_oster_delta.csv")

cat("\n--- S2: Oster delta ---\n")
print(s2_tbl)

latex_out <- c(latex_out,
  to_latex(s2_tbl,
           caption = "S2: Oster (2019) $\\delta$. Computed via OLS approximation on pooled data. $\\delta > 1$ indicates robustness by convention. Reference category: $\\gamma$ (Bridging). $R^2_{\\max} = \\min(1.3 \\times R^2_{\\text{full}}, 1)$.",
           label   = "tab:SI_S2_oster"))

# ==============================================================================
# S3 — CROSS-LAGGED (clase_{t-1} → trust_t | trust_{t-1})
# ==============================================================================
lg("\n--- S3: Cross-lagged specification ---")
#
# clase at wave t-1 predicts trust at wave t, controlling for lagged trust.
# Stronger evidence of temporal precedence than contemporaneous models.
# Uses only waves 2 and 3 (N ≈ 2×1297 = 2594 person-wave obs).

dt_lag <- dt %>%
  arrange(idencuesta, ola) %>%
  group_by(idencuesta) %>%
  mutate(
    clase_lag    = lag(clase),
    trust_lag    = lag(trust),
    trust_nh_lag = lag(trust_nh)
  ) %>%
  ungroup() %>%
  filter(ola > 1)   # drop wave 1 (no lag available)

lg("  Cross-lagged sample: ", nrow(dt_lag), " obs | ",
   length(unique(dt_lag$idencuesta)), " ids")

fit_crosslag <- function(outcome, lag_var, ref_clase, label) {

  d <- dt_lag %>%
    filter(!is.na(clase_lag), !is.na(.data[[outcome]]), !is.na(.data[[lag_var]])) %>%
    mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola),
      clase_f      = relevel(as.factor(clase_lag), ref = as.character(ref_clase))
    )

  ctrl <- ctrl_use[ctrl_use %in% names(d) &
                   sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  d <- drop_na(d, all_of(intersect(
    c(outcome, lag_var, "clase_f", "ola_f", ctrl), names(d)
  )))

  rhs <- paste(c("clase_f", lag_var, ctrl, "ola_f", "(1 | idencuesta_f)"),
               collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  lg("  Fitting cross-lag ", label, " | N=", nrow(d))
  m <- tryCatch(
    glmer(fml, data = d, family = binomial("probit"), nAGQ = NAGQ,
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 5e5))),
    error = function(e) {
      lg("  Cross-lag nAGQ=", NAGQ, " failed, retrying Laplace")
      glmer(fml, data = d, family = binomial("probit"), nAGQ = 1L,
            control = glmerControl(optimizer = "bobyqa"))
    }
  )

  avg_slopes(m, variables = "clase_f", type = "response") %>%
    as_tibble() %>%
    mutate(spec = label, outcome = outcome, model = "Cross-lagged",
           n_obs = nrow(d))
}

s3_r1 <- fit_crosslag("trust",    "trust_lag",    ref_clase = 2L,
                      label = "S3_trust_ref_gamma")
s3_r2 <- fit_crosslag("trust_nh", "trust_nh_lag", ref_clase = 2L,
                      label = "S3_nh_ref_gamma")

s3_ames <- bind_rows(s3_r1, s3_r2) %>%
  filter(!grepl("^ola|_lag", contrast, ignore.case = TRUE)) %>%
  label_contrasts()

write_csv(s3_ames, here::here("output", "SI_S3_crosslagged.csv"))
lg("Saved: SI_S3_crosslagged.csv")

cat("\n--- S3: Cross-lagged AMEs (clase_{t-1} → trust_t | trust_{t-1}) ---\n")
s3_ames %>%
  select(spec, outcome, contrast_label, estimate, std.error, p.value, sig) %>%
  print()

latex_out <- c(latex_out,
  to_latex(s3_ames %>% select(spec, outcome, contrast_label, estimate, std.error, p.value, sig),
           caption = "S3: Cross-lagged specification. Portfolio position at wave $t-1$ predicts trust at wave $t$, controlling for lagged trust. Waves 2 and 3 only. Reference: $\\gamma$ (Bridging).",
           label   = "tab:SI_S3_crosslag"))

# ==============================================================================
# S4 — ENTROPY-WEIGHTED (weights = p_max)
# ==============================================================================
lg("\n--- S4: Entropy-weighted regression ---")
#
# Each person-wave observation is weighted by its maximum posterior probability
# p_max = max(p_alpha, p_beta, p_gamma). Observations with uncertain
# classification contribute less to the likelihood.
# Weights are normalized to sum to N (preserves effective sample size).

fit_wtd <- function(outcome, ref_clase, label) {
  # Pooled probit with importance weights (p_max).
  # glm(family=binomial, weights=wt) treats wt as frequency weights,
  # which is the correct interpretation for importance-/entropy-weights
  # in a pooled model (unlike glmer where weights = number of Bernoulli trials).
  # Cluster-robust SEs (by idencuesta) via sandwich::vcovCL.

  d <- dt %>%
    filter(!is.na(clase), !is.na(.data[[outcome]]), !is.na(p_max)) %>%
    mutate(
      ola_f   = as.factor(ola),
      clase_f = relevel(as.factor(clase), ref = as.character(ref_clase)),
      wt      = p_max / mean(p_max, na.rm = TRUE)  # normalize to mean 1
    )

  ctrl <- ctrl_use[ctrl_use %in% names(d) &
                   sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  d <- drop_na(d, all_of(intersect(
    c(outcome, "clase_f", "ola_f", ctrl), names(d)
  )))

  rhs <- paste(c("clase_f", ctrl, "ola_f"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  lg("  Fitting weighted pooled probit ", label, " | N=", nrow(d),
     " | mean(p_max)=", round(mean(d$p_max), 3))

  m <- tryCatch(
    glm(fml, data = d, family = binomial("probit"), weights = wt),
    error = function(e) { lg("  Weighted glm failed: ", e$message); NULL }
  )
  if (is.null(m)) return(NULL)

  # Cluster-robust variance-covariance (by person)
  vcv <- tryCatch(
    sandwich::vcovCL(m, cluster = ~idencuesta),
    error = function(e) { lg("  vcovCL failed: ", e$message); NULL }
  )

  avg_slopes(m, variables = "clase_f", type = "response",
             vcov = if (!is.null(vcv)) vcv else TRUE) %>%
    as_tibble() %>%
    mutate(spec = label, outcome = outcome, model = "Entropy-weighted (pooled probit)",
           n_obs = nrow(d), mean_pmax = round(mean(d$p_max), 3))
}

s4_r1 <- fit_wtd("trust",    ref_clase = 2L, label = "S4_trust_ref_gamma")
s4_r2 <- fit_wtd("trust_nh", ref_clase = 2L, label = "S4_nh_ref_gamma")

s4_ames <- bind_rows(s4_r1, s4_r2) %>%
  filter(!grepl("^ola", contrast, ignore.case = TRUE)) %>%
  label_contrasts()

write_csv(s4_ames, here::here("output", "SI_S4_entropy_weighted.csv"))
lg("Saved: SI_S4_entropy_weighted.csv")

cat("\n--- S4: Entropy-weighted AMEs ---\n")
s4_ames %>%
  select(spec, outcome, contrast_label, estimate, std.error, p.value, sig) %>%
  print()

latex_out <- c(latex_out,
  to_latex(s4_ames %>% select(spec, outcome, contrast_label, estimate, std.error, p.value, sig),
           caption = "S4: Entropy-weighted pooled probit. Observations weighted by maximum posterior probability $p_{\\max}$ (normalized to mean 1). Cluster-robust standard errors by respondent. Reference: $\\gamma$ (Bridging).",
           label   = "tab:SI_S4_entropy"))

# ==============================================================================
# S5 — SOFT ASSIGNMENT (p_gamma, p_beta as continuous predictors)
# ==============================================================================
lg("\n--- S5: Soft assignment ---")
#
# Instead of modal classification, uses posterior probabilities directly as
# continuous predictors [0, 1]. Reference = p_alpha (implicit, omitted).
# The AME of p_gamma answers: "a 1-unit increase in the posterior probability
# of γ predicts how many pp change in trust?"
# This avoids the classification problem entirely.

fit_soft <- function(outcome, label) {

  d <- dt %>%
    filter(!is.na(p_gamma), !is.na(p_beta), !is.na(.data[[outcome]])) %>%
    mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola)
    )

  ctrl <- ctrl_use[ctrl_use %in% names(d) &
                   sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  d <- drop_na(d, all_of(intersect(
    c(outcome, "p_gamma", "p_beta", "ola_f", ctrl), names(d)
  )))

  rhs <- paste(c("p_gamma", "p_beta", ctrl, "ola_f", "(1 | idencuesta_f)"),
               collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  lg("  Fitting soft ", label, " | N=", nrow(d))
  m <- tryCatch(
    glmer(fml, data = d, family = binomial("probit"), nAGQ = NAGQ,
          control = glmerControl(optimizer = "bobyqa",
                                 optCtrl   = list(maxfun = 5e5))),
    error = function(e) {
      lg("  Soft nAGQ=", NAGQ, " failed, retrying Laplace")
      glmer(fml, data = d, family = binomial("probit"), nAGQ = 1L,
            control = glmerControl(optimizer = "bobyqa"))
    }
  )

  # AMEs for p_gamma and p_beta (the key predictors)
  ame <- avg_slopes(m, variables = c("p_gamma", "p_beta"), type = "response") %>%
    as_tibble() %>%
    mutate(
      spec    = label,
      outcome = outcome,
      model   = "Soft assignment",
      n_obs   = nrow(d),
      sig     = case_when(
        p.value < .001 ~ "***", p.value < .01 ~ "**",
        p.value < .05  ~ "*",   p.value < .10  ~ ".",
        TRUE           ~ ""
      )
    )

  list(model = m, ame = ame)
}

s5_r1 <- fit_soft("trust",    label = "S5_trust_soft")
s5_r2 <- fit_soft("trust_nh", label = "S5_nh_soft")

s5_ames <- bind_rows(s5_r1$ame, s5_r2$ame)

write_csv(s5_ames, here::here("output", "SI_S5_soft_assignment.csv"))
lg("Saved: SI_S5_soft_assignment.csv")

cat("\n--- S5: Soft assignment AMEs ---\n")
s5_ames %>%
  select(spec, outcome, term, estimate, std.error, p.value, sig) %>%
  print()

latex_out <- c(latex_out,
  to_latex(s5_ames %>% select(spec, outcome, term, estimate, std.error, p.value, sig),
           caption = "S5: Soft assignment specification. Posterior probabilities $p_{\\gamma}$ and $p_{\\beta}$ enter as continuous predictors (reference: $p_{\\alpha}$, omitted). AME reports marginal effect of a 1-unit increase in each posterior probability.",
           label   = "tab:SI_S5_soft"))

# ==============================================================================
# COMPARISON TABLE — AMEs across all specifications
# ==============================================================================
lg("\n--- Comparison table: all sensitivity specs ---")

# Standardize columns for comparison
std_cols <- function(df, spec_id) {
  df %>%
    mutate(sensitivity = spec_id) %>%
    select(sensitivity, outcome,
           contrast   = if ("contrast_label" %in% names(df)) "contrast_label" else "term",
           estimate, std.error, p.value, sig) %>%
    filter(!grepl("^ola", contrast, ignore.case = TRUE))
}

comparison <- bind_rows(
  std_cols(s1_ames %>% mutate(contrast_label = contrast_label), "S1_CRE"),
  std_cols(s3_ames, "S3_CrossLag"),
  std_cols(s4_ames, "S4_EntropyWtd"),
  s5_ames %>% mutate(sensitivity = "S5_Soft",
                     contrast = term) %>%
    select(sensitivity, outcome, contrast, estimate, std.error, p.value, sig)
)

write_csv(comparison, here::here("output", "SI_sensitivity_comparison.csv"))
lg("Saved: SI_sensitivity_comparison.csv")

cat("\n=== SENSITIVITY COMPARISON (p < .10) ===\n")
comparison %>%
  filter(p.value < .10) %>%
  mutate(est_sig = paste0(round(estimate, 3), sig)) %>%
  select(sensitivity, outcome, contrast, est_sig, std.error) %>%
  print(n = 40)

# ==============================================================================
# WRITE LaTeX TABLES
# ==============================================================================
latex_out <- c(
  "% === SI Sensitivity Tables — auto-generated by 03a_trust_sensitivity.R ===",
  "% Paste each table into your SI .tex file.\n",
  latex_out
)
writeLines(latex_out, here::here("output", "SI_sensitivity_tables.txt"))
lg("Saved: SI_sensitivity_tables.txt")

# ==============================================================================
# DONE
# ==============================================================================
lg("\n[03a_trust_sensitivity.R] DONE.")
lg("Outputs:")
lg("  SI_S1_CRE_mundlak.csv        — CRE/Mundlak AMEs + Wald test of mean terms")
lg("  SI_S2_oster_delta.csv        — Oster delta bounds")
lg("  SI_S3_crosslagged.csv        — Cross-lagged AMEs (clase_{t-1} -> trust_t)")
lg("  SI_S4_entropy_weighted.csv   — Entropy-weighted AMEs")
lg("  SI_S5_soft_assignment.csv    — Soft (posterior probability) AMEs")
lg("  SI_sensitivity_comparison.csv — All specs side by side")
lg("  SI_sensitivity_tables.txt    — LaTeX tables ready for SI")
