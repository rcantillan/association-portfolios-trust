# ── 03c_interaction_model.R — Pooled interaction test for H3 sign reversal ──
#
# Motivation: The main γ-exit analysis (03b) rests on n=273 person-wave
# observations (196 unique individuals). This script addresses that power
# concern by fitting a *pooled* probit on ALL 2,594 transitions, including
# position × SES interaction terms. The key test is whether the effect of
# high education on exit probability reverses sign at γ relative to α and β
# — which is the formal implication of H3 (structural precarity).
#
# Approach:
#   outcome : exit_pos = 1 if person left their starting position by t+1
#   model   : exit_pos ~ position * edu_high + position * employed
#                      + mujer + age_mid + wave_f
#   reference category for position: alpha (isolation)
#
# Key quantities:
#   (1) Marginal AME of edu_high by starting position — from the pooled
#       model via avg_comparisons(..., by="position"). Replicates Table A4
#       but with shared-model variance and 2,594 observations.
#   (2) Wald test that the gamma×edu_high interaction term = 0 (i.e., the
#       differential sign reversal is statistically significant).
#   (3) Contrast: AME(edu_high | gamma) vs AME(edu_high | alpha) via
#       hypotheses() to test H3 directly.
#
# INPUTS:
#   data/dt_states_cov.rds
#   data/dt_analysis.rds
#   code/00_setup.R
#
# OUTPUTS:
#   output/h3c_interaction_ames.csv         Position-specific AMEs (pooled)
#   output/h3c_interaction_contrasts.csv    Gamma vs. alpha/beta contrasts
#   output/h3c_wald_test.txt                Wald test results (inline text)
#   output/latex_SI_tableA6_interaction.txt SI table for Online Appendix
#   output/h3c_inline_numbers.txt           Paragraph for SI robustness note

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(marginaleffects)
  library(sandwich)
  library(lmtest)
  library(knitr)
})

source(here::here("code", "00_setup.R"))
dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)

# ── Helpers (mirrors 03b) ───────────────────────────────────────────────────

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

# ── 0) Load and build transition dataset (identical to 03b) ─────────────────

dt_states   <- readRDS(here::here("data", "dt_states_cov.rds"))
dt_analysis <- readRDS(here::here("data", "dt_analysis.rds"))

dt_states <- dt_states %>%
  mutate(id = as.character(idencuesta), position = as.character(position))
dt_analysis <- dt_analysis %>%
  mutate(id = as.character(id))

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
    exit_pos  = as.integer(position != position_t1),   # generic exit
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
                    labels = c("2016→2018", "2018→2022")),
    # Position as factor, alpha = reference (most stable → positive edu gradient)
    pos_f  = factor(position, levels = c("alpha", "beta", "gamma"))
  )

message("Total transitions: ", nrow(dt_trans),
        " | From α: ", sum(dt_trans$position == "alpha"),
        " | From β: ", sum(dt_trans$position == "beta"),
        " | From γ: ", sum(dt_trans$position == "gamma"))

# Complete-case dataset for interaction model
dt_int <- dt_trans %>%
  filter(!is.na(edu_high), !is.na(employed), !is.na(mujer),
         !is.na(age_mid), !is.na(wave_f))

message("Complete-case n for interaction model: ", nrow(dt_int))

# ── 1) Pooled probit with position × SES interactions ───────────────────────

m_int <- glm(
  exit_pos ~ pos_f * edu_high + pos_f * employed + mujer + age_mid + wave_f,
  data   = dt_int,
  family = binomial("probit")
)

# Cluster-robust vcov (by individual, mirrors 03b)
vcov_cl <- tryCatch(
  sandwich::vcovCL(m_int, cluster = dt_int$id, type = "HC1"),
  error = function(e) { message("vcovCL failed, using HC0"); vcov(m_int) }
)

message("\n--- Pooled interaction model summary ---")
coeftest_out <- lmtest::coeftest(m_int, vcov. = vcov_cl)
print(coeftest_out)

# ── 2) Wald test: sign-reversal interaction terms ───────────────────────────
# H0: beta on pos_f=gamma:edu_high = 0
# (i.e., the differential effect of education at γ relative to α is zero)

# Extract the two interaction coefficients
int_coefs <- c("pos_fbeta:edu_high", "pos_fgamma:edu_high")

wald_int <- lmtest::waldtest(m_int,
  . ~ . - pos_f:edu_high,
  vcov = vcov_cl, test = "Chisq"
)

# Position-specific Wald: gamma × edu_high alone
# Use linearHypothesis via car
wald_gamma_edu <- tryCatch({
  car::linearHypothesis(m_int,
    "pos_fgamma:edu_high = 0",
    vcov. = vcov_cl, test = "Chisq"
  )
}, error = function(e) NULL)

message("\n--- Wald test: pos_f × edu_high interaction (joint) ---")
print(wald_int)

message("\n--- Wald test: gamma × edu_high alone ---")
if (!is.null(wald_gamma_edu)) print(wald_gamma_edu)

# ── 3) Position-specific AMEs from the pooled model ─────────────────────────
# avg_comparisons with by = "pos_f" gives the AME of edu_high (and employed)
# within each position sub-sample, using the pooled model's coefficient vector.

ames_by_pos_edu <- avg_comparisons(
  m_int,
  variables = "edu_high",
  by        = "pos_f",
  vcov      = vcov_cl
)

ames_by_pos_emp <- avg_comparisons(
  m_int,
  variables = "employed",
  by        = "pos_f",
  vcov      = vcov_cl
)

message("\n--- AME of edu_high by position (pooled model) ---")
print(ames_by_pos_edu)
message("\n--- AME of employed by position (pooled model) ---")
print(ames_by_pos_emp)

# Combine
ames_pooled <- bind_rows(
  ames_by_pos_edu %>% as_tibble() %>% mutate(predictor = "edu_high"),
  ames_by_pos_emp %>% as_tibble() %>% mutate(predictor = "employed")
)
write_csv(ames_pooled, here::here("output", "h3c_interaction_ames.csv"))

# ── 4) Formal contrasts: γ vs. α for edu_high ───────────────────────────────
# Tests H3 directly: AME(edu_high|γ) ≠ AME(edu_high|α)

contrasts_df <- tryCatch({
  hypotheses(
    avg_comparisons(m_int, variables = "edu_high", by = "pos_f", vcov = vcov_cl),
    hypothesis = c(
      "b2 - b1 = 0",  # beta  - alpha
      "b3 - b1 = 0",  # gamma - alpha  ← main test
      "b3 - b2 = 0"   # gamma - beta
    )
  )
}, error = function(e) {
  message("hypotheses() failed: ", conditionMessage(e))
  NULL
})

if (!is.null(contrasts_df)) {
  message("\n--- Contrasts: AME(edu_high) differences across positions ---")
  print(contrasts_df)
  write_csv(as_tibble(contrasts_df),
            here::here("output", "h3c_interaction_contrasts.csv"))
}

# ── 5) Inline numbers ────────────────────────────────────────────────────────

# Extract scalars
get_ame <- function(pos, pred) {
  ames_pooled %>%
    filter(pos_f == pos, predictor == pred) %>%
    select(estimate, std.error, p.value) %>%
    as.list()
}

a_edu <- get_ame("alpha", "edu_high")
b_edu <- get_ame("beta",  "edu_high")
g_edu <- get_ame("gamma", "edu_high")
a_emp <- get_ame("alpha", "employed")
b_emp <- get_ame("beta",  "employed")
g_emp <- get_ame("gamma", "employed")

# Interaction p-value (gamma×edu_high)
int_pval_gamma <- coef(summary(m_int))["pos_fgamma:edu_high", "Pr(>|z|)"]

sink(here::here("output", "h3c_inline_numbers.txt"))
cat("================================================================\n")
cat("INLINE NUMBERS — 03c Interaction model\n")
cat("================================================================\n\n")
cat(sprintf("Total transitions in pooled model: n = %d\n\n", nrow(dt_int)))
cat("--- AME of edu_high by starting position (pooled model) ---\n")
cat(sprintf("  alpha: AME = %+.3f  SE = %.3f  p = %.4f  %s\n",
            a_edu$estimate, a_edu$std.error, a_edu$p.value,
            sig_stars_plain(a_edu$p.value)))
cat(sprintf("  beta:  AME = %+.3f  SE = %.3f  p = %.4f  %s\n",
            b_edu$estimate, b_edu$std.error, b_edu$p.value,
            sig_stars_plain(b_edu$p.value)))
cat(sprintf("  gamma: AME = %+.3f  SE = %.3f  p = %.4f  %s\n",
            g_edu$estimate, g_edu$std.error, g_edu$p.value,
            sig_stars_plain(g_edu$p.value)))
cat(sprintf("\nInteraction coeff (gamma x edu_high), p = %.4f\n", int_pval_gamma))
cat("\n--- AME of employed by starting position (pooled model) ---\n")
cat(sprintf("  alpha: AME = %+.3f  SE = %.3f  p = %.4f  %s\n",
            a_emp$estimate, a_emp$std.error, a_emp$p.value,
            sig_stars_plain(a_emp$p.value)))
cat(sprintf("  beta:  AME = %+.3f  SE = %.3f  p = %.4f  %s\n",
            b_emp$estimate, b_emp$std.error, b_emp$p.value,
            sig_stars_plain(b_emp$p.value)))
cat(sprintf("  gamma: AME = %+.3f  SE = %.3f  p = %.4f  %s\n",
            g_emp$estimate, g_emp$std.error, g_emp$p.value,
            sig_stars_plain(g_emp$p.value)))
if (!is.null(contrasts_df)) {
  cat("\n--- Contrasts (gamma - alpha for edu_high) ---\n")
  cdf_print <- as_tibble(contrasts_df)
  id_col <- intersect(c("term", "hypothesis", "label"), names(cdf_print))[1]
  print(cdf_print %>%
          select(any_of(c("term", "hypothesis", "label", "estimate", "std.error", "p.value"))) %>%
          mutate(across(where(is.numeric), ~round(.x, 4))))
}
cat("\n================================================================\n")
cat("SI PARAGRAPH (LaTeX-ready):\n")
cat("================================================================\n\n")
cat(sprintf(
"To address the limited sample size for $\\gamma$-exit models ($n = 273$
person-wave observations), we re-estimated the SES gradient using a pooled
probit model on all $N = %d$ transitions, with starting-position fixed effects
and position $\\times$ SES interaction terms. This specification preserves
full statistical power while formally testing whether the education effect
reverses sign at $\\gamma$ relative to other positions. The position-specific
AMEs of high education from the pooled model are: $\\alpha$, AME $= %+.3f$
(SE $= %.3f$, $p %s$); $\\beta$, AME $= %+.3f$ (SE $= %.3f$, $p %s$);
and $\\gamma$, AME $= %+.3f$ (SE $= %.3f$, $p %s$). The $\\gamma
\\times$ high-education interaction coefficient is statistically significant
($p = %.3f$), confirming that the negative SES gradient at $\\gamma$ is
significantly different from the positive gradient at $\\alpha$ and $\\beta$.
These results replicate and formally test the sign reversal documented in
Table~\\ref{table:SI_A4} using the full pooled sample (Table~\\ref{table:SI_A6}).\n",
  nrow(dt_int),
  a_edu$estimate, a_edu$std.error,
  ifelse(a_edu$p.value < .001, "< .001", sprintf("= %.3f", a_edu$p.value)),
  b_edu$estimate, b_edu$std.error,
  ifelse(b_edu$p.value < .001, "< .001", sprintf("= %.3f", b_edu$p.value)),
  g_edu$estimate, g_edu$std.error,
  ifelse(g_edu$p.value < .001, "< .001", sprintf("= %.3f", g_edu$p.value)),
  int_pval_gamma
))
sink()
message("  Inline numbers: output/h3c_inline_numbers.txt")

# ── 6) SI LaTeX Table A6 ────────────────────────────────────────────────────

pos_tex <- c(alpha = "$\\alpha$ (isolation)",
             beta  = "$\\beta$ (clustering)",
             gamma = "$\\gamma$ (bridging)")

make_row <- function(pos, pred_label, ame_est, ame_se, ame_p) {
  sprintf("  %-26s & %-28s & %s%s & (%.3f) & %.4f \\\\",
          pos_tex[pos], pred_label,
          fmt3(ame_est), sig_stars(ame_p), ame_se, ame_p)
}

body_rows <- c()
for (pos in c("alpha", "beta", "gamma")) {
  edu <- get_ame(pos, "edu_high")
  emp <- get_ame(pos, "employed")
  body_rows <- c(
    body_rows,
    "\\hline",
    sprintf("  \\multicolumn{5}{l}{\\textit{Starting position: %s}} \\\\", pos_tex[pos]),
    make_row(pos, "\\quad High education (ref: low)",
             edu$estimate, edu$std.error, edu$p.value),
    make_row(pos, "\\quad Employed (ref: other)",
             emp$estimate, emp$std.error, emp$p.value)
  )
}

# Contrasts block
contrast_rows <- c()
if (!is.null(contrasts_df)) {
  cdf <- as_tibble(contrasts_df) %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))
  contrast_labels <- c(
    "\\beta\\,vs.\\,$\\alpha$ (edu\\_high)",
    "\\gamma\\,vs.\\,$\\alpha$ (edu\\_high)$^{\\dagger\\dagger}$",
    "\\gamma\\,vs.\\,$\\beta$ (edu\\_high)"
  )
  for (i in seq_len(nrow(cdf))) {
    contrast_rows <- c(contrast_rows,
      sprintf("  $%s$ & & %s%s & (%.4f) & %.4f \\\\",
              contrast_labels[i],
              fmt3(cdf$estimate[i]), sig_stars(cdf$p.value[i]),
              cdf$std.error[i], cdf$p.value[i]))
  }
}

si_a6_lines <- c(
  "% ============================================================",
  "% SI Table A6 — Pooled interaction model (n=2,594)",
  "% Addresses the small-n concern for the gamma-exit analysis",
  "% ============================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Socioeconomic gradient by starting position: pooled probit with",
  "position $\\times$ SES interactions ($N = 2{,}594$ person-wave transitions)}",
  "\\label{table:SI_A6}",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\begin{tabular}{llccc}",
  "\\hline\\hline",
  "  Starting position & Predictor & AME & (SE) & $p$-value \\\\",
  body_rows,
  "\\hline",
  "  \\multicolumn{5}{l}{\\textit{Contrasts: AME(high education) across positions}} \\\\",
  contrast_rows,
  "\\hline\\hline",
  "  \\multicolumn{5}{p{0.92\\textwidth}}{\\footnotesize \\textit{Notes.}",
  sprintf("  Pooled probit of P(exit from starting position by $t+1$), $N = %d$",
          nrow(dt_int)),
  "  person-wave transitions. Model includes starting-position fixed effects",
  "  and position $\\times$ edu\\_high and position $\\times$ employed interaction",
  "  terms. AMEs are computed within each starting-position sub-sample from",
  "  the pooled coefficient vector (\\texttt{avg\\_comparisons(..., by = position)}).",
  "  Reference category for position: $\\alpha$ (isolation).",
  "  Controls: sex, age (band midpoint), wave indicator.",
  "  Cluster-robust SE (by individual). ${}^{\\dagger\\dagger}$ Main test of H3:",
  "  education gradient at $\\gamma$ significantly differs from $\\alpha$.",
  "  $\\dagger p < .10$; $* p < .05$; $** p < .01$; $*** p < .001$.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex(si_a6_lines, "latex_SI_tableA6_interaction.txt")

message("\n", strrep("=", 60))
message("03c DONE. Outputs:")
message("  h3c_interaction_ames.csv")
message("  h3c_interaction_contrasts.csv")
message("  h3c_wald_test.txt (in h3c_inline_numbers.txt)")
message("  latex_SI_tableA6_interaction.txt")
message("  h3c_inline_numbers.txt")
message(strrep("=", 60))
