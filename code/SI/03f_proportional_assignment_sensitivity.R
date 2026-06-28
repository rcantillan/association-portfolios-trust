# ── 03f_proportional_assignment_sensitivity.R ────────────────────────────────
#
# Sensitivity S8: Proportional (soft) class assignment
# -------------------------------------------------------
# The main models (03_trust_models.R) use MODAL assignment: each individual is
# hard-assigned to the class with the highest posterior probability. This is the
# standard "three-step" approach in LCA/LMM.
#
# A limitation is that modal assignment ignores classification uncertainty:
# an individual with posteriors (0.45, 0.44, 0.11) is treated identically to
# one with (0.99, 0.01, 0.00), although the former is much more ambiguous.
#
# This script implements PROPORTIONAL assignment (Yamaguchi, 2015; analogous to
# the "soft" three-step): instead of a hard 0/1 indicator for class membership,
# each person-wave observation contributes to all three class-outcome comparisons
# weighted by its posterior probability P(state | Y_1:T).
#
# Method:
#   1. Stack three copies of the data (one per class: alpha, beta, gamma).
#   2. Assign clase = k in copy k; set observation weight = p_k (posterior).
#   3. Fit the same RE probit as in 03_trust_models.R, using glmer(weights=w).
#   4. Compute AMEs as usual.
#
# Note on weights in glmer(): for binary outcomes with family=binomial(probit),
# the `weights` argument specifies the number of Bernoulli trials. Fractional
# weights between 0 and 1 are treated as fractional trials, which is numerically
# equivalent to importance weighting in this context and is how proportional
# three-step is typically implemented in software (e.g., Mplus).
#
# CONTROLS: Following the main specification, sociodemographic variables (age,
# sex, education) are NOT included. Controls are limited to swb, employed, couple
# (identical to CONTROLS_TRUST in 03_trust_models.R).
#
# Outputs:
#   output/tableA19_proportional_assignment.csv
#   output/tableA19_proportional_ame_comparison.csv
#
# References:
#   Bakk, Tekle, & Vermunt (2013). Estimating the association between LC
#     membership and external variables using bias-adjusted three-step LCA.
#   Yamaguchi (2015). Propensity score analysis using latent class analysis.
#
# ─────────────────────────────────────────────────────────────────────────────
source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(lme4); library(marginaleffects)
  library(sandwich); library(broom.mixed); library(stringr)
})
if (!requireNamespace("merDeriv", quietly = TRUE))
  install.packages("merDeriv", quiet = TRUE)
suppressPackageStartupMessages(library(merDeriv))

options(marginaleffects_safe = FALSE)

NAGQ <- 1L    # Laplace (nAGQ=1): sufficient for SI robustness check.
              # Main models use nAGQ=12; difference is negligible for AMEs.

# ── 1. Load data ───────────────────────────────────────────────────────────────
stop_if_missing(c(
  here::here("data", "dt_states_cov.rds"),
  here::here("data", "dt_analysis.rds")
))

dt_states <- readRDS(here::here("data", "dt_states_cov.rds")) %>% as_tibble()
dt_anal   <- readRDS(here::here("data", "dt_analysis.rds"))   %>% as_tibble()

# Merge trust and controls (same join as 03_trust_models.R)
dt_anal_join <- dt_anal %>%
  dplyr::rename(idencuesta = id) %>%
  dplyr::select(
    idencuesta, ola,
    trust, trust_nh,
    swb, employed, couple
  )

dt_states    <- dt_states    %>% dplyr::mutate(idencuesta = as.character(idencuesta))
dt_anal_join <- dt_anal_join %>% dplyr::mutate(idencuesta = as.character(idencuesta))

cols_new <- setdiff(names(dt_anal_join),
                    c("idencuesta", "ola", setdiff(names(dt_states),
                                                   c("idencuesta", "ola"))))
dt <- dplyr::left_join(
  dt_states,
  dt_anal_join %>% dplyr::select(dplyr::all_of(c("idencuesta", "ola", cols_new))),
  by = c("idencuesta", "ola")
)

# Verify posteriors are present
stopifnot(all(c("p_alpha", "p_beta", "p_gamma") %in% names(dt)))
message("Posteriors available. N=", nrow(dt), " obs | ",
        length(unique(dt$idencuesta)), " individuals")

# Sanity check: posteriors sum to ~1
post_sum <- dt %>%
  dplyr::filter(!is.na(p_alpha)) %>%
  dplyr::mutate(s = p_alpha + p_beta + p_gamma) %>%
  dplyr::summarise(min = min(s), max = max(s), mean = mean(s))
message("Posterior sums: min=", round(post_sum$min, 4),
        " max=", round(post_sum$max, 4),
        " mean=", round(post_sum$mean, 4))
if (post_sum$max > 1.01 | post_sum$min < 0.98)
  warning("Posterior probabilities do not sum to 1. Check dt_states_cov.rds.")

# Controls: same as main models — NO sociodemographics
CONTROLS_TRUST <- c("swb", "employed", "couple")
ctrl_use <- CONTROLS_TRUST[
  CONTROLS_TRUST %in% names(dt) &
  sapply(CONTROLS_TRUST, function(v) sum(!is.na(dt[[v]])) > 0)
]
message("Controls: ", paste(ctrl_use, collapse = ", "))

# ── 2. Build stacked proportional-assignment dataset ──────────────────────────
#
# For each person-wave observation, create 3 rows:
#   Row with clase=3 (alpha/Isolate),  weight = p_alpha
#   Row with clase=1 (beta/Closed),    weight = p_beta
#   Row with clase=2 (gamma/Bridging), weight = p_gamma
#
# clase coding matches canonical STATE_MAP (same as 03_trust_models.R):
#   alpha = Isolate   = clase 3
#   beta  = Closed    = clase 1   [reference when ref_clase=1]
#   gamma = Bridging  = clase 2   [reference when ref_clase=2]

dt_prop <- dplyr::bind_rows(
  dt %>% dplyr::mutate(clase = 3L, w_prop = p_alpha),   # alpha
  dt %>% dplyr::mutate(clase = 1L, w_prop = p_beta),    # beta
  dt %>% dplyr::mutate(clase = 2L, w_prop = p_gamma)    # gamma
) %>%
  dplyr::filter(!is.na(w_prop), w_prop > 1e-6)  # drop negligible-weight rows

message("Stacked proportional dataset: ", nrow(dt_prop), " rows | ",
        "effective N ≈ ", round(sum(dt_prop$w_prop), 0), " (should ≈ original N)")

# ── 3-5. Fit models and accumulate results directly ──────────────────────────
#
# Run 8 models inline (no helper function) to avoid S3 class issues.
# Each model: prep data → glmer → avg_slopes → extract key scalars → store row.

model_specs <- list(
  list(outcome="trust",    ref=2L, ctrl=FALSE, label="Prop_M1"),
  list(outcome="trust",    ref=2L, ctrl=TRUE,  label="Prop_M2"),
  list(outcome="trust_nh", ref=2L, ctrl=FALSE, label="Prop_M3"),
  list(outcome="trust_nh", ref=2L, ctrl=TRUE,  label="Prop_M4"),
  list(outcome="trust",    ref=1L, ctrl=FALSE, label="Prop_M5"),
  list(outcome="trust",    ref=1L, ctrl=TRUE,  label="Prop_M6"),
  list(outcome="trust_nh", ref=1L, ctrl=FALSE, label="Prop_M7"),
  list(outcome="trust_nh", ref=1L, ctrl=TRUE,  label="Prop_M8")
)

# Label map for contrast strings
label_map <- c(
  "1 - 2"="β vs γ", "3 - 2"="α vs γ",
  "2 - 1"="γ vs β", "3 - 1"="α vs β",
  "1 - 3"="β vs α", "2 - 3"="γ vs α"
)

results_rows <- list()
message("\n=== TABLE A19 — Proportional assignment sensitivity ===")

for (spec in model_specs) {
  message("\n  [PROP] Fitting: ", spec$label,
          " | outcome=", spec$outcome, " | ref=clase", spec$ref,
          " | controls=", spec$ctrl, " | nAGQ=", NAGQ)

  outcome <- spec$outcome
  ref_c   <- spec$ref
  lbl     <- spec$label
  ctrl    <- if (spec$ctrl) ctrl_use else character(0)

  d <- dt_prop
  d <- d[!is.na(d[[outcome]]) & !is.na(d$w_prop), ]
  d$idencuesta_f <- as.factor(d$idencuesta)
  d$ola_f        <- as.factor(d$ola)
  d$clase_f      <- relevel(as.factor(d$clase), ref = as.character(ref_c))

  # Drop NAs in required columns
  req_cols <- c(outcome, "clase_f", "ola_f", "w_prop", ctrl)
  req_cols <- intersect(req_cols, names(d))
  d <- d[complete.cases(d[, req_cols]), ]

  message("    N rows=", nrow(d),
          " | eff_N≈", round(sum(d$w_prop), 0),
          " | ids=", length(unique(d$idencuesta_f)))

  rhs <- paste(c("clase_f", ctrl, "ola_f", "(1 | idencuesta_f)"), collapse=" + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  m <- tryCatch(
    lme4::glmer(fml, data=d, family=binomial(link="probit"),
                weights=d$w_prop, nAGQ=NAGQ,
                control=lme4::glmerControl(optimizer="bobyqa",
                                           optCtrl=list(maxfun=5e5))),
    error = function(e) {
      message("  nAGQ=", NAGQ, " failed → Laplace")
      lme4::glmer(fml, data=d, family=binomial(link="probit"),
                  weights=d$w_prop, nAGQ=1L,
                  control=lme4::glmerControl(optimizer="bobyqa",
                                             optCtrl=list(maxfun=5e5)))
    }
  )

  sigma_u <- round(sqrt(as.numeric(lme4::VarCorr(m)$idencuesta_f)), 4)
  eff_n   <- round(sum(d$w_prop), 0)
  message("    sigma_u=", sigma_u)

  # Extract AMEs using avg_slopes; pull scalars directly to avoid S3 issues
  ame_obj <- marginaleffects::avg_slopes(m, variables="clase_f",
                                         type="response", vcov=TRUE)

  # Pull as named numeric vectors (safe, no S3 dependencies)
  est  <- as.numeric(ame_obj[["estimate"]])
  se   <- as.numeric(ame_obj[["std.error"]])
  pval <- as.numeric(ame_obj[["p.value"]])
  clo  <- as.numeric(ame_obj[["conf.low"]])
  chi  <- as.numeric(ame_obj[["conf.high"]])
  ctrs <- vapply(seq_len(nrow(ame_obj)),
                 function(i) as.character(ame_obj[["contrast"]][i]),
                 character(1))

  for (i in seq_along(est)) {
    ctr_str <- ctrs[i]
    ctr_lbl <- if (ctr_str %in% names(label_map)) label_map[[ctr_str]] else ctr_str
    results_rows[[length(results_rows) + 1L]] <- list(
      spec           = lbl,
      outcome        = outcome,
      ref_clase      = ref_c,
      with_controls  = spec$ctrl,
      contrast       = ctr_str,
      contrast_label = ctr_lbl,
      estimate       = est[i],
      std.error      = se[i],
      p.value        = pval[i],
      conf.low       = clo[i],
      conf.high      = chi[i],
      eff_n          = eff_n,
      sigma_u        = sigma_u
    )
  }
}

# Build table from list of lists — guaranteed plain data.frame
table_prop <- data.frame(
  spec           = vapply(results_rows, `[[`, character(1), "spec"),
  outcome        = vapply(results_rows, `[[`, character(1), "outcome"),
  ref_clase      = vapply(results_rows, `[[`, integer(1),  "ref_clase"),
  with_controls  = vapply(results_rows, `[[`, logical(1),  "with_controls"),
  contrast       = vapply(results_rows, `[[`, character(1), "contrast"),
  contrast_label = vapply(results_rows, `[[`, character(1), "contrast_label"),
  estimate       = vapply(results_rows, `[[`, numeric(1),  "estimate"),
  std.error      = vapply(results_rows, `[[`, numeric(1),  "std.error"),
  p.value        = vapply(results_rows, `[[`, numeric(1),  "p.value"),
  conf.low       = vapply(results_rows, `[[`, numeric(1),  "conf.low"),
  conf.high      = vapply(results_rows, `[[`, numeric(1),  "conf.high"),
  eff_n          = vapply(results_rows, `[[`, numeric(1),  "eff_n"),
  sigma_u        = vapply(results_rows, `[[`, numeric(1),  "sigma_u"),
  stringsAsFactors = FALSE
)

write.csv(table_prop, here::here("output", "tableA19_proportional_assignment.csv"),
          row.names = FALSE)
message("Saved: output/tableA19_proportional_assignment.csv")

# Display
cat("\n--- Table A19 (proportional assignment) ---\n")
table_prop$sig <- ifelse(table_prop$p.value<.001,"***",
                  ifelse(table_prop$p.value<.01,"**",
                  ifelse(table_prop$p.value<.05,"*",
                  ifelse(table_prop$p.value<.10,".","  "))))
print(table_prop[, c("spec","outcome","contrast_label","estimate","std.error","p.value","sig","eff_n","sigma_u")],
      row.names = FALSE, digits = 4)

# ── 6. Side-by-side comparison: modal vs proportional ────────────────────────
modal_path <- here::here("output", "ame_all_models.csv")
if (file.exists(modal_path)) {
  modal_raw <- read.csv(modal_path, stringsAsFactors = FALSE)
  modal <- modal_raw[
    modal_raw$table %in% c("Table4_ref_Bridging","Table5_ref_Closed") &
    !grepl("^ola|^Ola", modal_raw$contrast_label), ]
  modal <- modal[, c("spec","outcome","contrast_label","estimate","std.error","p.value","n_obs")]
  names(modal) <- c("spec_modal","outcome","contrast_label","est_modal","se_modal","p_modal","n_modal")
  modal$sig_modal <- ifelse(modal$p_modal<.001,"***",
                     ifelse(modal$p_modal<.01,"**",
                     ifelse(modal$p_modal<.05,"*",
                     ifelse(modal$p_modal<.10,".","  "))))

  # Map proportional spec labels to modal ones
  spec_map <- c(Prop_M1="T4_M1",Prop_M2="T4_M2",Prop_M3="T4_M3",Prop_M4="T4_M4",
                Prop_M5="T5_M5",Prop_M6="T5_M6",Prop_M7="T5_M7",Prop_M8="T5_M8")
  table_prop$spec_modal <- spec_map[table_prop$spec]
  table_prop$sig_prop   <- ifelse(table_prop$p.value<.001,"***",
                            ifelse(table_prop$p.value<.01,"**",
                            ifelse(table_prop$p.value<.05,"*",
                            ifelse(table_prop$p.value<.10,".","  "))))

  prop_sub <- table_prop[, c("spec","spec_modal","outcome","contrast_label",
                               "estimate","std.error","p.value","eff_n","sig_prop")]
  names(prop_sub)[names(prop_sub) %in% c("estimate","std.error","p.value","eff_n")] <-
    c("est_prop","se_prop","p_prop","n_prop")

  comparison <- merge(prop_sub, modal,
                      by.x = c("spec_modal","outcome","contrast_label"),
                      by.y = c("spec_modal","outcome","contrast_label"),
                      all.x = TRUE)
  comparison$delta     <- round(comparison$est_prop - comparison$est_modal, 4)
  comparison$pct_delta <- round((comparison$est_prop - comparison$est_modal) /
                                  abs(comparison$est_modal) * 100, 1)
  comparison$consistent <- sign(comparison$est_prop) == sign(comparison$est_modal)

  write.csv(comparison,
            here::here("output", "tableA19_proportional_ame_comparison.csv"),
            row.names = FALSE)
  message("Saved: output/tableA19_proportional_ame_comparison.csv")

  cat("\n--- Modal vs Proportional AME comparison ---\n")
  cat("(delta = prop - modal; pct_delta = % change relative to modal)\n\n")
  print(comparison[, c("spec","outcome","contrast_label",
                        "est_modal","sig_modal","est_prop","sig_prop",
                        "delta","pct_delta","consistent")],
        row.names = FALSE, digits = 4)

  n_consistent <- sum(comparison$consistent, na.rm = TRUE)
  n_total      <- sum(!is.na(comparison$consistent))
  cat("\nSign consistency:", n_consistent, "/", n_total, "contrasts\n")
  if (n_consistent == n_total)
    cat("  All contrasts have the same sign under both assignment rules.\n")
  else
    cat("  Some contrasts change sign — inspect manually.\n")

} else {
  message("  ame_all_models.csv not found. Run 03_trust_models.R first for comparison.")
}

message("\n[03f_proportional_assignment_sensitivity.R] done.")
message("  Output: output/tableA19_proportional_assignment.csv")
message("  Output: output/tableA19_proportional_ame_comparison.csv")
message("")
message("  Interpretation:")
message("  - Stable AMEs (similar magnitude, same sign, same significance)")
message("    across modal and proportional assignment rules indicate that")
message("    the portfolio-trust association is not an artifact of hard class")
message("    assignment and is robust to classification uncertainty.")
message("  - Proportional assignment does NOT resolve demographic confounding;")
message("    it addresses a different source of sensitivity (assignment rule).")
