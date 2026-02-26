# ==============================================================================
# 03_trust_models.R — Trust models + AMEs + γ−β contrast (DROP employment)
# ------------------------------------------------------------------------------
# Outputs:
#   output/trust_ame_table.csv
#   output/gamma_minus_beta_contrast.csv
#   output/model_sample_diagnostics.csv
#   output/trust_models_glmer_summaries.txt
# ------------------------------------------------------------------------------
# Key features (improved):
#   - robust input validation (required columns + binary outcomes)
#   - consistent sample construction per outcome (complete cases on model vars)
#   - guards against 1-level factors after filtering (position/ola/woman/education)
#   - cluster-robust SEs by id (CR2) to match Stata vce(cluster id)
#   - writes diagnostics + model summaries to disk for reproducibility
# ==============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(lme4)
  library(marginaleffects)
  library(readr)
  library(clubSandwich)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 0) Load data
# ------------------------------------------------------------------------------
stop_if_missing(c(here::here("data", "dt_states.rds")))
dt <- readRDS(here::here("data", "dt_states.rds")) |> as_tibble()

dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1) Validate required columns
# ------------------------------------------------------------------------------
required_cols <- c("id", "ola", "position", "trust", "trust_nh", "edad", "woman", "education")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in dt_states: ", paste(missing_cols, collapse = ", "),
    "\nAvailable columns: ", paste(names(dt), collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 2) Standardize types (avoid silent issues)
# ------------------------------------------------------------------------------
dt <- dt %>%
  mutate(
    id = as.factor(id),
    ola = as.factor(ola),
    # enforce α baseline; if position is numeric-state, coerce carefully:
    position = as.character(position),
    position = factor(position, levels = c("alpha", "beta", "gamma"))
  )

# Controls used in Stata MINUS employment
controls <- c("edad", "woman", "education")

# ------------------------------------------------------------------------------
# 3) Helper: enforce binary outcomes (0/1)
# ------------------------------------------------------------------------------
coerce_binary01 <- function(x, nm) {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x)) {
    lev <- levels(x)
    if (all(lev %in% c("0", "1"))) return(as.integer(as.character(x)))
    stop("Outcome `", nm, "` is a factor with levels: ", paste(lev, collapse = ", "),
         ". Convert to 0/1 before modeling.")
  }
  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (all(ux %in% c(0, 1))) return(as.integer(x))
    stop("Outcome `", nm, "` is not binary 0/1. Unique values (first 15): ",
         paste(head(ux, 15), collapse = ", "),
         if (length(ux) > 15) " ...")
  }
  stop("Outcome `", nm, "` has unsupported type: ", class(x)[1])
}

# ------------------------------------------------------------------------------
# 4) Cluster vcov (matches Stata vce(cluster id))
# ------------------------------------------------------------------------------
Vcov_cluster_id <- function(model) {
  clubSandwich::vcovCR(model, cluster = model@frame$id, type = "CR2")
}

# ------------------------------------------------------------------------------
# 5) Build model sample per outcome + diagnostics
# ------------------------------------------------------------------------------
build_model_sample <- function(outcome) {
  fml <- as.formula(
    paste0(outcome, " ~ position + ", paste(controls, collapse = " + "), " + ola + (1|id)")
  )
  vars <- all.vars(fml)

  d <- dt %>%
    select(all_of(vars)) %>%
    mutate(
      # coerce binary outcome (throws if not 0/1)
      "{outcome}" := coerce_binary01(.data[[outcome]], outcome),
      # ensure factors (and drop unused levels after filtering)
      position = droplevels(position),
      ola = droplevels(ola)
    ) %>%
    filter(complete.cases(.)) %>%
    mutate(
      position = droplevels(position),
      ola = droplevels(ola),
      # woman/education: keep as-is; if they are categorical, convert to factor
      woman = if (n_distinct(woman) <= 10) woman else woman,
      education = if (n_distinct(education) <= 20) education else education
    )

  # convert common categorical controls to factors (safe)
  # (Stata uses i.woman i.education, so we factorize if they have few unique values)
  if ("woman" %in% names(d) && n_distinct(d$woman) <= 5) d$woman <- as.factor(d$woman)
  if ("education" %in% names(d) && n_distinct(d$education) <= 50) d$education <- as.factor(d$education)

  # post-factor drop levels
  d$woman <- if (is.factor(d$woman)) droplevels(d$woman) else d$woman
  d$education <- if (is.factor(d$education)) droplevels(d$education) else d$education

  # Guardrails: need variation where it matters
  if (nrow(d) == 0) stop("No complete cases for outcome: ", outcome)
  if (nlevels(d$id) < 2) stop("After filtering, `id` has <2 levels for: ", outcome)
  if (nlevels(d$position) < 2) stop("After filtering, `position` has <2 levels for: ", outcome,
                                   " (check dt_states construction).")
  if (nlevels(d$ola) < 2) stop("After filtering, `ola` has <2 levels for: ", outcome)

  # If factors woman/education collapsed to 1 level, drop them (otherwise contrasts error)
  dropped <- character(0)
  rhs <- c("position", controls, "ola")

  if (is.factor(d$woman) && nlevels(d$woman) < 2) { rhs <- setdiff(rhs, "woman"); dropped <- c(dropped, "woman") }
  if (is.factor(d$education) && nlevels(d$education) < 2) { rhs <- setdiff(rhs, "education"); dropped <- c(dropped, "education") }

  fml2 <- as.formula(paste0(outcome, " ~ ", paste(rhs, collapse = " + "), " + (1|id)"))

  list(
    data = d,
    fml = fml2,
    dropped = dropped,
    diag = tibble(
      outcome = outcome,
      n_rows = nrow(d),
      n_ids = nlevels(d$id),
      pos_levels = nlevels(d$position),
      ola_levels = nlevels(d$ola),
      woman_levels = if (is.factor(d$woman)) nlevels(d$woman) else n_distinct(d$woman),
      edu_levels = if (is.factor(d$education)) nlevels(d$education) else n_distinct(d$education),
      dropped_terms = if (length(dropped) == 0) "" else paste(dropped, collapse = ",")
    )
  )
}

# ------------------------------------------------------------------------------
# 6) Fit models
# ------------------------------------------------------------------------------
fit_re_probit <- function(outcome) {
  ms <- build_model_sample(outcome)

  if (nchar(ms$diag$dropped_terms) > 0) {
    message("Dropping 1-level covariates for ", outcome, ": ", ms$diag$dropped_terms)
  }

  model <- glmer(
    ms$fml, data = ms$data,
    family = binomial(link = "probit"),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )

  attr(model, "model_diag") <- ms$diag
  model
}

m_g <- fit_re_probit("trust")
m_n <- fit_re_probit("trust_nh")

diag_tbl <- bind_rows(attr(m_g, "model_diag"), attr(m_n, "model_diag"))
write_csv(diag_tbl, here::here("output", "model_sample_diagnostics.csv"))
print(diag_tbl)

# Save model summaries (helps debugging / exact replication)
summ_txt <- c(
  "=== Model: trust (generalized) ===",
  capture.output(summary(m_g)),
  "",
  "=== Model: trust_nh (neighborhood) ===",
  capture.output(summary(m_n))
)
writeLines(summ_txt, here::here("output", "trust_models_glmer_summaries.txt"))

# ------------------------------------------------------------------------------
# 7) AMEs (probability scale) + γ−β contrast
# ------------------------------------------------------------------------------
ame_g <- avg_slopes(m_g, variables = "position", type = "response", vcov = Vcov_cluster_id) %>%
  as_tibble() %>% mutate(outcome = "generalized_trust")

ame_n <- avg_slopes(m_n, variables = "position", type = "response", vcov = Vcov_cluster_id) %>%
  as_tibble() %>% mutate(outcome = "neighborhood_trust")

ame_tbl <- bind_rows(ame_g, ame_n)
write_csv(ame_tbl, here::here("output", "trust_ame_table.csv"))

cmp_g <- avg_comparisons(
  m_g,
  variables = "position",
  comparisons = list(gamma_minus_beta = c("gamma", "beta")),
  type = "response",
  vcov = Vcov_cluster_id
) %>% as_tibble() %>% mutate(outcome = "generalized_trust")

cmp_n <- avg_comparisons(
  m_n,
  variables = "position",
  comparisons = list(gamma_minus_beta = c("gamma", "beta")),
  type = "response",
  vcov = Vcov_cluster_id
) %>% as_tibble() %>% mutate(outcome = "neighborhood_trust")

cmp_tbl <- bind_rows(cmp_g, cmp_n)
write_csv(cmp_tbl, here::here("output", "gamma_minus_beta_contrast.csv"))

print(ame_tbl)
print(cmp_tbl)

# ------------------------------------------------------------------------------
# 8) Optional console diagnostics (quick parity check)
# ------------------------------------------------------------------------------
diag_levels <- function(outcome) {
  ms <- build_model_sample(outcome)
  d <- ms$data

  cat("\n=== Outcome:", outcome, "===\n")
  cat("N rows:", nrow(d), " | N ids:", nlevels(d$id), "\n")
  cat("position levels:", paste(levels(d$position), collapse = ", "), "\n")
  cat("ola levels:", paste(levels(d$ola), collapse = ", "), "\n")

  cat("\nCounts by position:\n")
  print(count(d, position, sort = TRUE))

  cat("\nCounts by ola:\n")
  print(count(d, ola, sort = TRUE))

  cat("\nCounts by position x ola:\n")
  print(count(d, position, ola) %>% arrange(ola, position))
}

diag_levels("trust")
diag_levels("trust_nh")