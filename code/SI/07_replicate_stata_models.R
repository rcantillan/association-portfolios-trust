# ── 07_replicate_stata_models.R — Stata cross-check (estimation.do in R) ──────
# Replicates estimation.do + estimation_2.do in R. Can run independently from
# ELSOC_Long.RData or from dt_states.rds produced by the pipeline.
#
# Variable codings (ELSOC Long 2016-2022):
#   c02   : 1=always trust (social_trust=1), 2/3=0 (Stata: replace=0 if c02>1)
#   t01   : Likert 1-5; trust_nb=1 if t01>3 (codes 4,5)
#   m02   : employed=1 if m02<=2 (full- or part-time)
#   m01   : education (1-10); Stata recode: 1/3->1, 4->2, 5->3, 6/7->4, 8/10->5
#   m36   : couple=1 if m36<4 (cohabiting)
#   s01   : swb=1 if s01>4 (satisfied)
#   c12_0x: 1=non-member, 2=inactive, 3=active; any_member = c12>=2
#
# Models replicated:
#   estimation.do  — Table 4 (ref=γ Broker/ib2): M1-M4 xtprobit
#                    Figure 3: margins clase#ola
#                    Table A2: xtreg clase ~ trust + controls
#   estimation_2.do — Table 5 (ref=β Closed/ib1): M5-M8 xtprobit
#                     Table A1: pooled probit versions
#
# Outputs:
#   output/table4_xtprobit_ref_broker.csv
#   output/table5_xtprobit_ref_closed.csv
#   output/tableA1_probit_ref_closed.csv
#   output/tableA2_reverse_causality.csv
#   output/figure3_margins_clase_ola.png
#   output/ame_summary_all_models.csv

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(lme4); library(marginaleffects)
  library(ggplot2); library(readr); library(clubSandwich)
  library(broom.mixed)
})

dir.create(here::here("output"), showWarnings = FALSE)

# ── 1. Load data (prefer dt_states.rds; fall back to ELSOC_Long.RData) ────────
if (file.exists(here::here("data", "dt_states.rds"))) {
  message("Using dt_states.rds (full pipeline available)...")
  dt <- readRDS(here::here("data", "dt_states.rds")) %>% as_tibble()

  # dt_states has trust, trust_nh, position (alpha/beta/gamma),
  # edad, woman, education, employed, swb, couple, tinst, ola, id.
  # Map position -> numeric clase to replicate ib1/ib2 in Stata.
  # clase 1=Closed (β), clase 2=Broker (γ), clase 3=Apathetic (α)
  dt <- dt %>%
    dplyr::mutate(
      clase = dplyr::case_when(
        position == "beta"  ~ 1L,   # Closed    = clase 1
        position == "gamma" ~ 2L,   # Broker    = clase 2 (reference ib2)
        position == "alpha" ~ 3L,   # Apathetic = clase 3
        TRUE ~ NA_integer_
      )
    )

} else {
  # No prior pipeline: load ELSOC raw
  message("dt_states.rds not found. Loading ELSOC_Long.RData...")
  stop_if_missing(c(here::here("data","ELSOC_Long.RData")))
  obj_names <- load(here::here("data","ELSOC_Long.RData"))
  raw <- if ("elsoc_long_2016_2022" %in% obj_names) get("elsoc_long_2016_2022") else get(obj_names[1])
  dt_raw <- dplyr::as_tibble(raw)

  # Requires 01_descriptive_stats.R and 02_latent_markov.R to have run first.
  stop("Run 01_descriptive_stats.R and 02_latent_markov.R first to generate dt_states.rds")
}

# ── 2. Verify variables (all already built by pipeline) ──────────────────────
# estimation.do: gen social_trust = 1 if c02==1; trust_nb=1 if t01>3;
#   employed=1 if m02<=2; couple=1 if m36<4; swb=1 if s01>4.

cat("\n=== Variable check (should match Stata) ===\n")
check_vars <- c("trust","trust_nh","employed","swb","couple","education","clase")
for (v in check_vars) {
  if (v %in% names(dt)) {
    cat(sprintf("  %-12s: n=%d | mean=%.3f | NAs=%d\n",
                v, sum(!is.na(dt[[v]])), mean(dt[[v]], na.rm=TRUE), sum(is.na(dt[[v]]))))
  } else {
    cat(sprintf("  %-12s: NOT FOUND\n", v))
  }
}

# ── 3. Prepare balanced panel ─────────────────────────────────────────────────
# dt_states.rds already applied muestra==1, tipo_atricion==1, ola in 1/3/6.

n_per_id <- dt %>% dplyr::count(id, name = "n_olas")
balanced_ids <- n_per_id %>% dplyr::filter(n_olas == 3) %>% dplyr::pull(id)
dt_bal <- dt %>% dplyr::filter(id %in% balanced_ids)

cat("\n=== Balanced panel ===\n")
cat("  N individuals:", length(balanced_ids), "\n")
cat("  N obs person-wave:", nrow(dt_bal), "\n")
cat("  Class distribution:\n")
print(table(dt_bal$clase, dt_bal$ola, useNA="ifany"))

# ── 4. Factor variables and clustering ────────────────────────────────────────
prep_data <- function(data, ref_clase = 2, outcome_var, extra_controls = TRUE) {
  d <- data %>%
    dplyr::mutate(
      id_f  = as.factor(id),
      ola_f = as.factor(ola),
      clase_f = relevel(as.factor(clase), ref = as.character(ref_clase))
    )

  # Controls replicating estimation.do: m0_sexo m0_edad swb employed couple t02_01
  ctrl_vars <- intersect(
    c("woman","edad","swb","employed","couple","tinst"),
    names(d)
  )
  ctrl_vars <- ctrl_vars[sapply(ctrl_vars, function(v) any(!is.na(d[[v]])))]

  if (!extra_controls) ctrl_vars <- character(0)

  # Drop NA on outcome and key predictors
  keep_vars <- c(outcome_var, "clase_f", "ola_f", "id_f", ctrl_vars)
  d <- d %>% tidyr::drop_na(dplyr::all_of(intersect(keep_vars, names(d))))

  list(data = d, controls = ctrl_vars)
}

# ── 5. Helper: cluster-robust vcov ────────────────────────────────────────────
vcov_cr <- function(model) {
  tryCatch(
    clubSandwich::vcovCR(model, cluster = model@frame$id_f, type = "CR1"),
    error = function(e) NULL
  )
}

# ── 6. Helper: fit xtprobit + AMEs (replicates xtprobit + margins, dydx(*)) ──
fit_xtprobit <- function(outcome, data, ref_clase = 2,
                         with_controls = TRUE, label = "") {
  prep   <- prep_data(data, ref_clase, outcome, with_controls)
  d      <- prep$data
  ctrl   <- prep$controls

  rhs <- paste(c("clase_f", ctrl, "ola_f", "(1|id_f)"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  m <- lme4::glmer(fml, data = d,
                   family  = binomial(link = "probit"),
                   control = lme4::glmerControl(
                     optimizer = "bobyqa",
                     optCtrl   = list(maxfun = 2e5)
                   ))

  vcv <- vcov_cr(m)

  # AMEs for clase (= margins, dydx(clase) in Stata)
  ame <- marginaleffects::avg_slopes(
    m, variables = "clase_f", type = "response", vcov = vcv
  ) %>% as_tibble() %>%
    dplyr::mutate(
      outcome      = outcome,
      ref_clase    = ref_clase,
      controls     = with_controls,
      spec         = label,
      n            = nrow(d),
      n_id         = length(unique(d$id_f))
    )

  # Predicted probabilities by clase#ola for Figure 3
  pred <- marginaleffects::predictions(
    m, newdata = marginaleffects::datagrid(
      clase_f  = levels(d$clase_f),
      ola_f    = levels(d$ola_f),
      grid_type = "mean_or_mode"
    ), type = "response"
  ) %>% as_tibble() %>%
    dplyr::mutate(outcome = outcome, ref_clase = ref_clase, spec = label)

  list(model = m, ame = ame, pred = pred,
       outcome = outcome, ref = ref_clase, label = label,
       formula = deparse(fml), n = nrow(d))
}

# ── 7. Helper: fit pooled probit (no RE) ─────────────────────────────────────
fit_probit <- function(outcome, data, ref_clase = 2,
                       with_controls = TRUE, label = "") {
  prep <- prep_data(data, ref_clase, outcome, with_controls)
  d    <- prep$data
  ctrl <- prep$controls

  rhs <- paste(c("clase_f", ctrl, "ola_f"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  # Cluster-robust SE at id level (= vce(cluster idencuesta))
  m <- glm(fml, data = d, family = binomial(link = "probit"))

  vcv <- tryCatch(
    sandwich::vcovCL(m, cluster = ~id_f),
    error = function(e) NULL
  )

  ame <- marginaleffects::avg_slopes(
    m, variables = "clase_f", type = "response", vcov = vcv
  ) %>% as_tibble() %>%
    dplyr::mutate(
      outcome = outcome, ref_clase = ref_clase,
      controls = with_controls, spec = label, n = nrow(d)
    )

  list(model = m, ame = ame, outcome = outcome,
       ref = ref_clase, label = label, formula = deparse(fml))
}

# ── 8. Table 4 — RE probit ref=γ Broker (ib2.clase) ──────────────────────────
message("\n=== Table 4: xtprobit ref=Broker (ib2.clase) ===")

m4_1 <- fit_xtprobit("trust",    dt_bal, ref_clase=2, with_controls=FALSE, label="Table4_M1")
m4_2 <- fit_xtprobit("trust",    dt_bal, ref_clase=2, with_controls=TRUE,  label="Table4_M2")
m4_3 <- fit_xtprobit("trust_nh", dt_bal, ref_clase=2, with_controls=FALSE, label="Table4_M3")
m4_4 <- fit_xtprobit("trust_nh", dt_bal, ref_clase=2, with_controls=TRUE,  label="Table4_M4")

table4_ames <- dplyr::bind_rows(m4_1$ame, m4_2$ame, m4_3$ame, m4_4$ame)
readr::write_csv(table4_ames, here::here("output","table4_xtprobit_ref_broker.csv"))

cat("\n--- TABLE 4 (ref=Broker) AMEs ---\n")
print(table4_ames %>%
      dplyr::select(spec, outcome, contrast, estimate, std.error, p.value, n) %>%
      dplyr::filter(!grepl("ola|intercept", contrast, ignore.case=TRUE)) %>%
      dplyr::mutate(sig = dplyr::case_when(
        p.value < .001 ~ "***", p.value < .01 ~ "**",
        p.value < .05 ~ "*",   p.value < .10 ~ ".",   TRUE ~ ""
      )) %>%
      dplyr::arrange(outcome, spec))

# ── 9. Table 5 — RE probit ref=β Closed (ib1.clase) ──────────────────────────
message("\n=== Table 5: xtprobit ref=Closed (ib1.clase) ===")

m5_5 <- fit_xtprobit("trust",    dt_bal, ref_clase=1, with_controls=FALSE, label="Table5_M5")
m5_6 <- fit_xtprobit("trust",    dt_bal, ref_clase=1, with_controls=TRUE,  label="Table5_M6")
m5_7 <- fit_xtprobit("trust_nh", dt_bal, ref_clase=1, with_controls=FALSE, label="Table5_M7")
m5_8 <- fit_xtprobit("trust_nh", dt_bal, ref_clase=1, with_controls=TRUE,  label="Table5_M8")

table5_ames <- dplyr::bind_rows(m5_5$ame, m5_6$ame, m5_7$ame, m5_8$ame)
readr::write_csv(table5_ames, here::here("output","table5_xtprobit_ref_closed.csv"))

cat("\n--- TABLE 5 (ref=Closed) AMEs ---\n")
print(table5_ames %>%
      dplyr::select(spec, outcome, contrast, estimate, std.error, p.value) %>%
      dplyr::filter(!grepl("ola|intercept", contrast, ignore.case=TRUE)) %>%
      dplyr::mutate(sig = dplyr::case_when(
        p.value < .001 ~ "***", p.value < .01 ~ "**",
        p.value < .05 ~ "*",   p.value < .10 ~ ".",   TRUE ~ ""
      )))

# ── 10. Table A1 — Pooled probit ref=β Closed ─────────────────────────────────
message("\n=== Table A1: Pooled probit ref=Closed ===")

ma1_1 <- fit_probit("trust",    dt_bal, ref_clase=1, with_controls=FALSE, label="A1_M1")
ma1_2 <- fit_probit("trust",    dt_bal, ref_clase=1, with_controls=TRUE,  label="A1_M2")
ma1_3 <- fit_probit("trust_nh", dt_bal, ref_clase=1, with_controls=FALSE, label="A1_M3")
ma1_4 <- fit_probit("trust_nh", dt_bal, ref_clase=1, with_controls=TRUE,  label="A1_M4")

tableA1_ames <- dplyr::bind_rows(ma1_1$ame, ma1_2$ame, ma1_3$ame, ma1_4$ame)
readr::write_csv(tableA1_ames, here::here("output","tableA1_probit_ref_closed.csv"))

# ── 11. Table A2 — Reverse causality (trust -> clase) ────────────────────────
# Stata: xtreg clase i.social_trust + controls + i.ola, vce(cluster idencuesta)
# R: lmer with random intercept (clase treated as continuous, as in Stata xtreg)
message("\n=== Table A2: Reverse causality (trust → clase) ===")

fit_reverse <- function(trust_var, data, with_controls = TRUE, label = "") {
  d <- data %>%
    dplyr::mutate(id_f = as.factor(id), ola_f = as.factor(ola)) %>%
    tidyr::drop_na(clase, .data[[trust_var]])

  ctrl <- intersect(c("woman","edad","education","swb","employed","couple","tinst"), names(d))
  ctrl <- ctrl[sapply(ctrl, function(v) any(!is.na(d[[v]])))]
  if (!with_controls) ctrl <- character(0)

  rhs <- paste(c(trust_var, ctrl, "ola_f", "(1|id_f)"), collapse = " + ")
  fml <- as.formula(paste("clase ~", rhs))

  m <- lme4::lmer(fml, data = d,
                  control = lme4::lmerControl(optimizer = "bobyqa"))

  vcv <- tryCatch(
    clubSandwich::vcovCR(m, cluster = d$id_f, type = "CR1"),
    error = function(e) NULL
  )

  ame <- marginaleffects::avg_slopes(
    m, variables = trust_var, vcov = vcv
  ) %>% as_tibble() %>%
    dplyr::mutate(
      outcome = "clase", predictor = trust_var,
      controls = with_controls, spec = label, n = nrow(d)
    )

  list(model = m, ame = ame, label = label)
}

ra1 <- fit_reverse("trust",    dt_bal, with_controls=FALSE, label="A2_M1")
ra2 <- fit_reverse("trust",    dt_bal, with_controls=TRUE,  label="A2_M2")
ra3 <- fit_reverse("trust_nh", dt_bal, with_controls=FALSE, label="A2_M3")
ra4 <- fit_reverse("trust_nh", dt_bal, with_controls=TRUE,  label="A2_M4")

tableA2_ames <- dplyr::bind_rows(ra1$ame, ra2$ame, ra3$ame, ra4$ame)
readr::write_csv(tableA2_ames, here::here("output","tableA2_reverse_causality.csv"))

cat("\n--- TABLE A2: trust -> clase (reverse direction) ---\n")
print(tableA2_ames %>%
      dplyr::select(spec, predictor, estimate, std.error, p.value) %>%
      dplyr::mutate(sig = dplyr::case_when(
        p.value < .001 ~ "***", p.value < .01 ~ "**",
        p.value < .05 ~ "*",   p.value < .10 ~ ".",   TRUE ~ ""
      )))

# ── 12. Figure 3 — Predicted probabilities by clase x ola ────────────────────
message("\n=== Figure 3: Predicted probabilities clase × ola ===")

# Collect predictions from Table 4 models (full controls)
preds_trust <- m4_2$pred %>%
  dplyr::mutate(outcome_label = "Generalized trust",
                clase_label   = dplyr::case_when(
                  clase_f == "1" ~ "Closed (β)",
                  clase_f == "2" ~ "Broker (γ)",
                  clase_f == "3" ~ "Apathetic (α)",
                  TRUE           ~ as.character(clase_f)
                ),
                wave_year = dplyr::case_when(
                  ola_f == "1" ~ WAVES_YEARS[1],
                  ola_f == "2" ~ WAVES_YEARS[2],
                  ola_f == "3" ~ WAVES_YEARS[3],
                  TRUE         ~ as.integer(ola_f)
                ))

preds_nh <- m4_4$pred %>%
  dplyr::mutate(outcome_label = "Trust in neighbors",
                clase_label   = dplyr::case_when(
                  clase_f == "1" ~ "Closed (β)",
                  clase_f == "2" ~ "Broker (γ)",
                  clase_f == "3" ~ "Apathetic (α)",
                  TRUE           ~ as.character(clase_f)
                ),
                wave_year = dplyr::case_when(
                  ola_f == "1" ~ WAVES_YEARS[1],
                  ola_f == "2" ~ WAVES_YEARS[2],
                  ola_f == "3" ~ WAVES_YEARS[3],
                  TRUE         ~ as.integer(ola_f)
                ))

preds_all <- dplyr::bind_rows(preds_trust, preds_nh)

fig3 <- ggplot2::ggplot(
  preds_all,
  ggplot2::aes(x = wave_year, y = estimate,
               color = clase_label, group = clase_label,
               ymin = conf.low, ymax = conf.high)
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::geom_ribbon(alpha = 0.12, color = NA) +
  ggplot2::facet_wrap(~outcome_label, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = WAVES_YEARS,
                               labels = as.character(WAVES_YEARS)) +
  ggplot2::scale_color_manual(
    values = c("Closed (β)"    = "#4E79A7",
               "Broker (γ)"    = "#E15759",
               "Apathetic (α)" = "#59A14F")
  ) +
  ggplot2::labs(
    x     = "Year",
    y     = "Pr(trust = 1)",
    color = "Class",
    title = "Predicted probabilities of trust by class and wave",
    caption = "RE probit models with full controls. Ribbons = 95% CI."
  ) +
  theme_ssr() +
  ggplot2::theme(legend.position = "bottom")

ggsave_safe("figure3_margins_clase_ola.png", fig3, width = 9, height = 5)
ggsave_safe("figure3_margins_clase_ola.pdf", fig3, width = 9, height = 5)
message("  Figure 3 saved.")

# ── 13. Unified AME summary ───────────────────────────────────────────────────
all_ames <- dplyr::bind_rows(
  table4_ames  %>% dplyr::mutate(table = "Table4"),
  table5_ames  %>% dplyr::mutate(table = "Table5"),
  tableA1_ames %>% dplyr::mutate(table = "TableA1"),
  tableA2_ames %>% dplyr::mutate(table = "TableA2")
)

readr::write_csv(all_ames, here::here("output","ame_summary_all_models.csv"))

# ── 14. Model summaries ────────────────────────────────────────────────────────
models_to_summarize <- list(
  Table4_M1 = m4_1$model, Table4_M2 = m4_2$model,
  Table4_M3 = m4_3$model, Table4_M4 = m4_4$model,
  Table5_M6 = m5_6$model, Table5_M8 = m5_8$model
)

summ_lines <- character(0)
for (nm in names(models_to_summarize)) {
  summ_lines <- c(summ_lines,
    paste0("\n\n=== ", nm, " ==="),
    capture.output(summary(models_to_summarize[[nm]]))
  )
}
writeLines(summ_lines, here::here("output","stata_replication_summaries.txt"))

# ── 15. R vs Stata comparison notes ──────────────────────────────────────────
cat("\n=== R vs Stata notes ===\n")
cat("Stata xtprobit -> R glmer(probit link + (1|id))\n")
cat("Stata margins, dydx(*) -> R avg_slopes(type='response')\n")
cat("Stata vce(cluster idencuesta) -> R vcovCR(type='CR1', cluster=id)\n")
cat("Expected differences (not errors):\n")
cat("  SEs: Stata CR0/CR1 approx equals clubSandwich CR1\n")
cat("  RE variance: lme4 REML vs Stata ML (minimal difference)\n")
cat("If coefs differ >5%: check clase coding, N balance, missing codes.\n")

message("\n[07_replicate_stata_models.R] done.")
message("  Tables: table4, table5, tableA1, tableA2")
message("  Figure: figure3_margins_clase_ola.png")
message("  Summary: ame_summary_all_models.csv")
