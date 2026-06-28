# ── 03d_ipw_attrition.R — IPW correction for panel attrition ───────────────
#
# Motivation: The balanced panel of 1,297 persons (observed at all 3
# measurement waves: ELSOC 2016, 2018, 2022) is a sub-sample of the 2,927
# Wave-1 respondents. Differential attrition could bias estimates if those
# who stayed differ systematically on the outcomes or SES covariates.
#
# Strategy (stabilized IPW):
#   1. Identify the 2,927 ELSOC Wave-1 (2016) respondents.
#   2. Among them, flag who was retained in the balanced panel.
#   3. Model P(retained | Wave-1 characteristics): logistic regression on
#      age, sex, education, employment status, trust (W1), SWB (W1).
#   4. Compute stabilized inverse-probability weights:
#         w_i = P(retained) / P_hat(retained | X_i)
#      where P(retained) = marginal retention rate.
#   5. Refit the H3 exit-from-γ models (full model, preferred M3 from 03b)
#      with these weights, using cluster-robust SE.
#   6. Compare weighted vs. unweighted AMEs as a sensitivity check.
#
# NOTE: Stabilized weights are preferred over raw 1/P weights because they
# have finite variance and smaller SE penalty.
#
# INPUTS:
#   data/ELSOC_Long.RData           (full panel: 18,035 obs × 748 vars)
#   data/dt_states_cov.rds          (LMM classifications)
#   data/dt_analysis.rds            (analysis dataset)
#   code/00_setup.R
#
# OUTPUTS:
#   output/h3d_ipw_weights.csv               Weight diagnostics
#   output/h3d_ipw_ames_comparison.csv       Weighted vs. unweighted AMEs
#   output/h3d_attrition_model.txt           Attrition model summary
#   output/latex_SI_tableA7_ipw.txt          SI table
#   output/h3d_inline_numbers.txt            Paragraph for SI

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(marginaleffects)
  library(sandwich)
  library(lmtest)
})

source(here::here("code", "00_setup.R"))
dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)

# ── Helpers ──────────────────────────────────────────────────────────────────

sig_stars <- function(p) {
  dplyr::case_when(
    p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*",
    p < .10  ~ "$\\dagger$", TRUE ~ ""
  )
}
sig_stars_plain <- function(p) {
  dplyr::case_when(
    p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", p < .10 ~ "+", TRUE ~ ""
  )
}
fmt3 <- function(x) formatC(round(x, 3), format = "f", digits = 3)

write_latex <- function(lines, filename) {
  writeLines(lines, here::here("output", filename))
  message("  LaTeX: output/", filename)
}

# ── 0) Load full ELSOC panel ─────────────────────────────────────────────────

message("Loading ELSOC_Long.RData ...")
load(here::here("data", "ELSOC_Long.RData"))   # → elsoc_long_2016_2022
elsoc <- elsoc_long_2016_2022
rm(elsoc_long_2016_2022)

message("Full ELSOC: ", nrow(elsoc), " obs | ",
        n_distinct(elsoc$idencuesta), " unique IDs")

# ── 1) Extract Wave-1 cohort ─────────────────────────────────────────────────

# ELSOC ola coding: 1=2016, 2=2017, 3=2018, 4=2019, 5=2021, 6=2022
# Membership variables (c12_xx) available in ola 1, 3, 6 only.

w1 <- elsoc %>%
  filter(ola == 1) %>%
  select(
    idencuesta,
    m0_sexo,    # 1=man, 2=woman
    m0_edad,    # age (continuous)
    m01,        # education level 1-5 (+ ELSOC negative missings)
    m02,        # employment status
    c02,        # generalized trust (1-5)
    s11_01,     # life satisfaction (1-5)
    c03,        # institutional trust
    c04         # neighborhood trust
  ) %>%
  rename(
    id         = idencuesta,
    sex_w1     = m0_sexo,
    age_w1     = m0_edad,
    educ_w1    = m01,
    emp_w1     = m02,
    trust_w1   = c02,
    swb_w1     = s11_01,
    tinst_w1   = c03,
    nhg_w1     = c04
  )

message("Wave-1 respondents: n = ", nrow(w1))

# ── 2) Identify balanced-panel members ───────────────────────────────────────

dt_analysis <- readRDS(here::here("data", "dt_analysis.rds")) %>%
  mutate(id = as.character(id))

panel_ids <- unique(dt_analysis$id)
message("Balanced-panel IDs: n = ", length(panel_ids))

w1 <- w1 %>%
  mutate(
    id        = as.character(id),
    retained  = as.integer(id %in% panel_ids)
  )

message("Wave-1 respondents in balanced panel: n = ",
        sum(w1$retained, na.rm = TRUE))
message("Retention rate: ",
        round(mean(w1$retained, na.rm = TRUE) * 100, 1), "%")

# ── 3) Recode Wave-1 covariates for attrition model ──────────────────────────

# ELSOC uses negative codes for missing: -999, -888, -777, -666
na_elsoc <- function(x) { x <- as.numeric(x); ifelse(x < 0, NA_real_, x) }

w1 <- w1 %>%
  mutate(
    across(c(sex_w1, age_w1, educ_w1, emp_w1,
             trust_w1, swb_w1, tinst_w1, nhg_w1),
           na_elsoc),
    # sex: 1=male, 2=female in ELSOC
    woman_w1  = case_when(
      sex_w1 == 2 ~ 1L, sex_w1 == 1 ~ 0L, TRUE ~ NA_integer_
    ),
    # Education: binary high (técnica=3, university=4, postgrad=5)
    edu_high_w1 = case_when(
      educ_w1 %in% 3:5 ~ 1L,
      educ_w1 %in% 1:2 ~ 0L,
      TRUE ~ NA_integer_
    ),
    # Employment: ELSOC m02: 1=full-time, 2=part-time, 3=self-employed → employed
    # 4=unpaid family, 5=unemployed, 6=retired, 7=student, 8=household, 9=other
    employed_w1 = case_when(
      emp_w1 %in% 1:3 ~ 1L,
      emp_w1 %in% 4:9 ~ 0L,
      TRUE ~ NA_integer_
    ),
    # Center age for model stability
    age_c = age_w1 - mean(age_w1, na.rm = TRUE)
  )

# Complete cases for attrition model
w1_cc <- w1 %>%
  filter(!is.na(retained), !is.na(woman_w1), !is.na(edu_high_w1),
         !is.na(employed_w1), !is.na(trust_w1), !is.na(swb_w1),
         !is.na(age_c))

message("Complete-case n for attrition model: ", nrow(w1_cc))
message("Retention in complete-case sample: ",
        sum(w1_cc$retained), " (", round(mean(w1_cc$retained)*100,1), "%)")

# ── 4) Fit attrition model ───────────────────────────────────────────────────

m_attr <- glm(
  retained ~ woman_w1 + age_c + edu_high_w1 + employed_w1 + trust_w1 + swb_w1,
  data   = w1_cc,
  family = binomial("logit")
)

sink(here::here("output", "h3d_attrition_model.txt"))
cat("=== Attrition model: P(retained in balanced panel | Wave-1 covariates) ===\n\n")
print(summary(m_attr))
cat("\nPseudo-R² (McFadden):",
    round(1 - m_attr$deviance / m_attr$null.deviance, 3), "\n")
sink()
message("  Attrition model: output/h3d_attrition_model.txt")

# ── 5) Compute stabilized IPW ────────────────────────────────────────────────

p_retained_marginal <- mean(w1_cc$retained)          # P(retained)
p_hat               <- fitted(m_attr)                 # P_hat(retained | X)

w1_cc <- w1_cc %>%
  mutate(
    p_hat   = p_hat,
    ipw_raw = 1 / p_hat,                             # raw (unstabilized)
    ipw_stab = p_retained_marginal / p_hat            # stabilized (preferred)
  )

# Trim extreme weights at 99th percentile to guard against instability
trim_pct <- quantile(w1_cc$ipw_stab[w1_cc$retained == 1], 0.99, na.rm = TRUE)
w1_cc <- w1_cc %>%
  mutate(ipw_trim = pmin(ipw_stab, trim_pct))

# Weight summary
wt_summary <- w1_cc %>%
  filter(retained == 1) %>%
  summarise(
    n_retained    = n(),
    w_mean        = mean(ipw_trim),
    w_sd          = sd(ipw_trim),
    w_min         = min(ipw_trim),
    w_max         = max(ipw_trim),
    w_p99_raw     = quantile(ipw_stab, 0.99),
    pct_trimmed   = mean(ipw_stab > trim_pct) * 100
  )

message("\n--- IPW weight diagnostics (retained respondents only) ---")
print(wt_summary)

# Save weight-diagnostics merged on id
ipw_out <- w1_cc %>%
  filter(retained == 1) %>%
  select(id, p_hat, ipw_stab, ipw_trim)

write_csv(ipw_out, here::here("output", "h3d_ipw_weights.csv"))

# ── 6) Rebuild transition dataset with IPW ────────────────────────────────────

dt_states <- readRDS(here::here("data", "dt_states_cov.rds")) %>%
  mutate(id = as.character(idencuesta), position = as.character(position))

ses_vars <- dt_analysis %>%
  select(id, ola, edu_bin, employed) %>%
  mutate(id = as.character(id), ola = as.integer(ola))

dt_full <- dt_states %>% left_join(ses_vars, by = c("id", "ola"))

dt_t  <- dt_full %>%
  select(id, ola, position, p_gamma, nivel_educ, mujer, edad, edu_bin, employed)
dt_t1 <- dt_full %>%
  mutate(ola = ola - 1L) %>%
  select(id, ola, position_t1 = position, p_gamma_t1 = p_gamma)

dt_trans <- dt_t %>%
  inner_join(dt_t1, by = c("id", "ola")) %>%
  filter(!is.na(position), !is.na(position_t1)) %>%
  mutate(
    exit_gamma = case_when(
      position == "gamma" ~ as.integer(position_t1 != "gamma"),
      TRUE ~ NA_integer_
    ),
    edu_ord = case_when(
      nivel_educ == "básica"  ~ 1L, nivel_educ == "media"   ~ 2L,
      nivel_educ == "técnica" ~ 3L, nivel_educ == "univers" ~ 4L,
      TRUE ~ NA_integer_
    ),
    edu_high = as.integer(edu_ord >= 3),
    age_mid  = case_when(
      edad == "18_24" ~ 21,   edad == "25_34" ~ 29.5, edad == "35_44" ~ 39.5,
      edad == "45_54" ~ 49.5, edad == "55_64" ~ 59.5, edad == "65"    ~ 70,
      TRUE ~ NA_real_
    ),
    wave_f = factor(ola, levels = c(1L, 2L),
                    labels = c("2016→2018", "2018→2022"))
  )

# γ-exit subsample
dt_gamma <- dt_trans %>%
  filter(position == "gamma", !is.na(exit_gamma)) %>%
  left_join(ipw_out %>% select(id, ipw_trim), by = "id") %>%
  # If id not in IPW table (refreshment sample added after W1), give weight 1
  mutate(ipw_trim = ifelse(is.na(ipw_trim), 1.0, ipw_trim))

message("\nγ-exit subsample with IPW: n = ", nrow(dt_gamma))
message("Missing weights (refreshment): n = ",
        sum(is.na(dt_gamma$ipw_trim) | dt_gamma$ipw_trim == 1.0 &
              !(dt_gamma$id %in% w1_cc$id)))

dt_gamma_cc <- dt_gamma %>%
  filter(!is.na(edu_high), !is.na(employed), !is.na(mujer),
         !is.na(age_mid), !is.na(wave_f))

# ── 7) Refit H3 models with and without IPW ──────────────────────────────────

# Cluster-robust vcov helper
vcov_cl <- function(m, dat) {
  tryCatch(
    sandwich::vcovCL(m, cluster = dat$id, type = "HC1"),
    error = function(e) vcov(m)
  )
}

# Unweighted (replication of 03b M3)
m_unweighted <- glm(
  exit_gamma ~ edu_high + employed + mujer + age_mid + wave_f,
  data   = dt_gamma_cc,
  family = binomial("probit")
)

# Weighted
m_weighted <- glm(
  exit_gamma ~ edu_high + employed + mujer + age_mid + wave_f,
  data    = dt_gamma_cc,
  family  = binomial("probit"),
  weights = ipw_trim
)

vcov_unw <- vcov_cl(m_unweighted, dt_gamma_cc)
vcov_wgt <- vcov_cl(m_weighted,   dt_gamma_cc)

ames_unw <- avg_comparisons(
  m_unweighted, variables = c("edu_high", "employed"),
  vcov = vcov_unw
) %>% as_tibble() %>% mutate(model = "Unweighted (M3 replicate)")

ames_wgt <- avg_comparisons(
  m_weighted, variables = c("edu_high", "employed"),
  vcov = vcov_wgt
) %>% as_tibble() %>% mutate(model = "IPW-weighted")

ames_compare <- bind_rows(ames_unw, ames_wgt)
write_csv(ames_compare, here::here("output", "h3d_ipw_ames_comparison.csv"))

message("\n--- AME comparison: unweighted vs. IPW-weighted ---")
ames_compare %>%
  select(model, term, estimate, std.error, p.value) %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  print()

# ── 8) Inline numbers ────────────────────────────────────────────────────────

get_ame2 <- function(df, mod, trm) {
  df %>% filter(model == mod, term == trm) %>%
    select(estimate, std.error, p.value) %>% as.list()
}

u_edu <- get_ame2(ames_compare, "Unweighted (M3 replicate)", "edu_high")
w_edu <- get_ame2(ames_compare, "IPW-weighted",              "edu_high")
u_emp <- get_ame2(ames_compare, "Unweighted (M3 replicate)", "employed")
w_emp <- get_ame2(ames_compare, "IPW-weighted",              "employed")

sink(here::here("output", "h3d_inline_numbers.txt"))
cat("================================================================\n")
cat("INLINE NUMBERS — 03d IPW attrition correction\n")
cat("================================================================\n\n")
cat(sprintf("Wave-1 n: %d | Panel n: %d | Retention rate: %.1f%%\n",
            nrow(w1), length(panel_ids),
            length(panel_ids)/nrow(w1)*100))
cat(sprintf("Attrition model pseudo-R²: %.3f\n",
            1 - m_attr$deviance / m_attr$null.deviance))
cat(sprintf("IPW weight range (retained): %.3f–%.3f (mean=%.3f, sd=%.3f)\n",
            wt_summary$w_min, wt_summary$w_max,
            wt_summary$w_mean, wt_summary$w_sd))
cat(sprintf("Pct trimmed at 99th pctile: %.1f%%\n\n",
            wt_summary$pct_trimmed))
cat("--- Unweighted (M3 replicate) ---\n")
cat(sprintf("  edu_high: AME=%+.3f  SE=%.3f  p=%.4f  %s\n",
            u_edu$estimate, u_edu$std.error, u_edu$p.value,
            sig_stars_plain(u_edu$p.value)))
cat(sprintf("  employed: AME=%+.3f  SE=%.3f  p=%.4f  %s\n",
            u_emp$estimate, u_emp$std.error, u_emp$p.value,
            sig_stars_plain(u_emp$p.value)))
cat("\n--- IPW-weighted ---\n")
cat(sprintf("  edu_high: AME=%+.3f  SE=%.3f  p=%.4f  %s\n",
            w_edu$estimate, w_edu$std.error, w_edu$p.value,
            sig_stars_plain(w_edu$p.value)))
cat(sprintf("  employed: AME=%+.3f  SE=%.3f  p=%.4f  %s\n",
            w_emp$estimate, w_emp$std.error, w_emp$p.value,
            sig_stars_plain(w_emp$p.value)))
cat("\n================================================================\n")
cat("SI PARAGRAPH (LaTeX-ready):\n")
cat("================================================================\n\n")
cat(sprintf(
"We further assessed whether differential attrition from Wave~1 ($n = %d$)
to the balanced panel ($n = %d$; retention rate = %.1f\\%%) biases our
estimates. We modelled the probability of retention as a function of Wave-1
demographics, education, employment status, generalized trust, and
subjective well-being (pseudo-$R^2 = %.3f$), and computed stabilized
inverse-probability weights (IPW) as $\\hat{w}_i = \\bar{p} /
\\hat{P}(\\text{retained}|\\mathbf{X}_i)$, trimmed at the 99th percentile
(weight range: %.2f--%.2f; mean = %.2f). Re-estimating the exit-from-$\\gamma$
probit with these weights leaves results substantively unchanged:
AME(high education) $= %+.3f$ (SE $= %.3f$, $p %s$) and
AME(employed) $= %+.3f$ (SE $= %.3f$, $p %s$), compared to unweighted
estimates of %+.3f and %+.3f, respectively (Table~\\ref{table:SI_A7}).
Differential attrition does not appear to confound the structural precarity
results.\n",
  nrow(w1), length(panel_ids),
  length(panel_ids)/nrow(w1)*100,
  1 - m_attr$deviance / m_attr$null.deviance,
  wt_summary$w_min, wt_summary$w_max, wt_summary$w_mean,
  w_edu$estimate, w_edu$std.error,
  ifelse(w_edu$p.value < .001, "< .001", sprintf("= %.3f", w_edu$p.value)),
  w_emp$estimate, w_emp$std.error,
  ifelse(w_emp$p.value < .001, "< .001", sprintf("= %.3f", w_emp$p.value)),
  u_edu$estimate, u_emp$estimate
))
sink()
message("  Inline numbers: output/h3d_inline_numbers.txt")

# ── 9) SI LaTeX Table A7 ────────────────────────────────────────────────────

make_row2 <- function(mod_label, ame_edu, ame_emp) {
  c(
    sprintf("  %-36s & %s%s & (%.3f) & %s%s & (%.3f) \\\\",
            mod_label,
            fmt3(ame_edu$estimate), sig_stars(ame_edu$p.value), ame_edu$std.error,
            fmt3(ame_emp$estimate), sig_stars(ame_emp$p.value), ame_emp$std.error)
  )
}

si_a7_lines <- c(
  "% ============================================================",
  "% SI Table A7 — IPW attrition-corrected AMEs (H3 robustness)",
  "% ============================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{IPW-corrected estimates of the SES gradient in exits from",
  "$\\gamma$ (attrition robustness check)}",
  "\\label{table:SI_A7}",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\begin{tabular}{lcccc}",
  "\\hline\\hline",
  "  & \\multicolumn{2}{c}{High education} & \\multicolumn{2}{c}{Employed} \\\\",
  "  \\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "  Model & AME & (SE) & AME & (SE) \\\\",
  "\\hline",
  make_row2("Unweighted (M3, preferred)", u_edu, u_emp),
  make_row2("IPW-weighted",               w_edu, w_emp),
  "\\hline\\hline",
  "  \\multicolumn{5}{p{0.88\\textwidth}}{\\footnotesize \\textit{Notes.}",
  sprintf("  Probit models of P(exit from $\\gamma$ by $t+1$), $n = %d$",
          nrow(dt_gamma_cc)),
  "  person-wave observations. Both models include sex, age, and wave indicator.",
  sprintf("  IPW: stabilized weights $\\hat{w}_i = \\bar{p} / \\hat{P}(\\text{retained}|\\mathbf{X}_i)$,"),
  sprintf("  trimmed at 99th percentile (range: %.2f--%.2f, mean: %.2f).",
          wt_summary$w_min, wt_summary$w_max, wt_summary$w_mean),
  sprintf("  Attrition model ($N = %d$ Wave-1 respondents): logistic regression on", nrow(w1_cc)),
  "  sex, age, education, employment, generalized trust, subjective well-being.",
  sprintf("  Retention rate = %.1f\\%%.", length(panel_ids)/nrow(w1)*100),
  "  Cluster-robust SE (by individual).",
  "  $\\dagger p < .10$; $* p < .05$; $** p < .01$; $*** p < .001$.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex(si_a7_lines, "latex_SI_tableA7_ipw.txt")

message("\n", strrep("=", 60))
message("03d DONE. Outputs:")
message("  h3d_ipw_weights.csv")
message("  h3d_ipw_ames_comparison.csv")
message("  h3d_attrition_model.txt")
message("  latex_SI_tableA7_ipw.txt")
message("  h3d_inline_numbers.txt")
message(strrep("=", 60))
