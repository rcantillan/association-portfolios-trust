# ==============================================================================
# replicate_estimation_do_UPDATED.R — Replicate estimation.do (Feb 28, 2024) in R
# ------------------------------------------------------------------------------
# INPUTS:
#   (preferred)   data/ELSOC_Long.dta    + data/data.long.dta
#   (fallback)    /mnt/data/ELSOC_Long.dta + /mnt/data/data.long.dta
#
# OUTPUTS (written to /output):
#   - output/table4_R.csv            (RE probit AMEs; 4 models)
#   - output/tableA1_R.csv           (Pooled probit AMEs; 4 models)
#   - output/tableA2_R.csv           (Reverse relation: trust -> class; 4 models)
#   - output/figure3_R.png / .pdf    (marginsplot-style predicted margins)
#   - output/diagnostics_estimation_do.txt
#   - output/balanced_panel_for_debug.csv
#
# Key updates vs prior:
#   (1) Rescale age (m0_edad) to centered decades -> reduces "nearly unidentifiable"
#   (2) Cluster vector for vcovCR comes from model@frame (no mismatch)
#   (3) Includes reverse relationship (Table A2) in R
#
# NOTE:
#   - Figure 3 is *predicted margins / predicted probabilities* (Stata: margins clase#ola + marginsplot)
#   - Tables are *average marginal effects* (Stata: margins, dydx(*))
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(haven)          # read_dta
  library(lme4)           # glmer / lmer
  library(sandwich)       # vcovCL
  library(clubSandwich)   # vcovCR
  library(marginaleffects)
})

DIR_OUT <- here::here("output")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 0) Locate inputs (data/ first, then /mnt/data fallback)
# ------------------------------------------------------------------------------
pick_path <- function(...) {
  p <- c(...)
  p <- p[file.exists(p)]
  if (length(p) == 0) stop("Input file not found in any expected location.")
  p[1]
}

PATH_ELSOC <- pick_path(
  here::here("data", "ELSOC_Long.dta"),
  "/mnt/data/ELSOC_Long.dta"
)
PATH_DL <- pick_path(
  here::here("data", "data.long.dta"),
  "/mnt/data/data.long.dta"
)

# ------------------------------------------------------------------------------
# 1) Helpers — replicate Stata mvdecode + Stata missing-comparison quirks
# ------------------------------------------------------------------------------
mvdecode_all <- function(df, codes = c(-888, -999, -666)) {
  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  df[num_cols] <- lapply(df[num_cols], function(x) {
    x[x %in% codes] <- NA
    x
  })
  df
}

# Stata:
# gen social_trust = 1 if c02 == 1
# replace social_trust = 0 if c02 > 1
# In Stata, missing(.) > 1 is TRUE -> missing becomes 0
code_social_trust <- function(c02) {
  c02 <- suppressWarnings(as.numeric(c02))
  ifelse(!is.na(c02) & c02 == 1, 1L, 0L)
}

# Stata:
# gen trust_nb = 1 if t01 > 3
# replace trust_nb = 0 if t01 <= 3
# In Stata, missing(.) > 3 is TRUE; missing(.) <= 3 is FALSE -> missing stays 1
code_trust_nb <- function(t01) {
  t01 <- suppressWarnings(as.numeric(t01))
  ifelse(is.na(t01) | t01 > 3, 1L, 0L)
}

# employed = 1 if m02 <= 2 ; else 0 ; missing -> 0 (because . > 2 is TRUE)
code_employed <- function(m02) {
  m02 <- suppressWarnings(as.numeric(m02))
  ifelse(!is.na(m02) & m02 <= 2, 1L, 0L)
}

# couple = 1 if m36 < 4 ; else 0 ; missing -> 0 (because . >= 4 is TRUE)
code_couple <- function(m36) {
  m36 <- suppressWarnings(as.numeric(m36))
  ifelse(!is.na(m36) & m36 < 4, 1L, 0L)
}

# swb = 1 if s01 > 4 ; else 0 ; missing -> 1 (because . > 4 is TRUE)
code_swb <- function(s01) {
  s01 <- suppressWarnings(as.numeric(s01))
  ifelse(is.na(s01) | s01 > 4, 1L, 0L)
}

# education recode (Stata):
# recode m01 (1/3=1) (4=2) (5=3) (6/7=4) (8/10=5) (.=.), gen(educ)
code_educ5 <- function(m01) {
  m01 <- suppressWarnings(as.numeric(m01))
  dplyr::case_when(
    is.na(m01)    ~ NA_integer_,
    m01 %in% 1:3  ~ 1L,
    m01 == 4      ~ 2L,
    m01 == 5      ~ 3L,
    m01 %in% 6:7  ~ 4L,
    m01 %in% 8:10 ~ 5L,
    TRUE          ~ NA_integer_
  )
}

stars <- function(p) {
  dplyr::case_when(p < .01 ~ "***", p < .05 ~ "**", p < .10 ~ "*", TRUE ~ "")
}

# Cluster-robust VCOV for glmer/lmer (cluster pulled from model frame)
get_vcov_mer <- function(model, cluster_var = "id_f") {
  cl <- tryCatch(model@frame[[cluster_var]], error = function(e) NULL)
  if (is.null(cl)) return(NULL)

  vcv <- tryCatch(
    clubSandwich::vcovCR(model, cluster = cl, type = "CR1"),
    error = function(e) NULL
  )
  if (is.null(vcv)) return(NULL)

  ev <- tryCatch(
    eigen(as.matrix(vcv), symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(ev) || any(ev < -1e-8)) return(NULL)
  vcv
}

# ------------------------------------------------------------------------------
# 2) Load + merge (Stata: merge 1:1 idencuesta ola using data.long.dta; keep _merge==3)
# ------------------------------------------------------------------------------
elsoc <- haven::read_dta(PATH_ELSOC)
dl    <- haven::read_dta(PATH_DL)

# Check uniqueness like merge 1:1
dup_elsoc <- elsoc %>% count(idencuesta, ola) %>% filter(n > 1)
dup_dl    <- dl    %>% count(idencuesta, ola) %>% filter(n > 1)
if (nrow(dup_elsoc) > 0) stop("ELSOC_Long.dta not unique on (idencuesta, ola).")
if (nrow(dup_dl) > 0)    stop("data.long.dta not unique on (idencuesta, ola).")

dt_merged <- inner_join(elsoc, dl, by = c("idencuesta", "ola"))

# mvdecode (-888, -999, -666)
dt_merged <- mvdecode_all(dt_merged, codes = c(-888, -999, -666))

# ------------------------------------------------------------------------------
# 3) Create variables exactly as Stata does (estimation.do)
# ------------------------------------------------------------------------------
need <- c("c02","t01","s01","s02","s03","m01","m02","m36",
          "muestra","tipo_atricion","clase","m0_sexo","m0_edad","t02_01","ola","idencuesta")
miss <- setdiff(need, names(dt_merged))
if (length(miss) > 0) stop("Missing required variables after merge: ", paste(miss, collapse=", "))

dt <- dt_merged %>%
  mutate(
    social_trust = code_social_trust(c02),
    trust_nb     = code_trust_nb(t01),
    life_aval    = (as.numeric(s01) + as.numeric(s02)) / 2,
    educ         = code_educ5(m01),
    unemployed   = ifelse(!is.na(m02) & as.numeric(m02) == 6, 1L, 0L),
    employed     = code_employed(m02),
    couple       = code_couple(m36),
    swb          = code_swb(s01),
    sh           = ifelse(!is.na(s03) & as.numeric(s03) >= 3, 1L, 0L)
  )

# ------------------------------------------------------------------------------
# 4) Filters (Stata):
#   keep if muestra==1
#   keep if tipo_atricion==1
#   keep if ola==1 | ola==3 | ola==6
# ------------------------------------------------------------------------------
dt <- dt %>%
  filter(muestra == 1, tipo_atricion == 1, ola %in% c(1, 3, 6))

# ------------------------------------------------------------------------------
# 5) Balance sample exactly like Stata:
#   probit social_trust ib2.clase i.m0_sexo m0_edad i.swb i.employed i.couple i.t02_01 i.ola, vce(cluster idencuesta)
#   gen in_model = e(sample); keep if in_model==1
#   keep if n_group==3
# ------------------------------------------------------------------------------
vars_esample <- c("social_trust","clase","m0_sexo","m0_edad","swb","employed","couple","t02_01","ola","idencuesta")
dt$in_model <- complete.cases(dt[, vars_esample])
dt <- dt %>% filter(in_model) %>% select(-in_model)

dt <- dt %>%
  group_by(idencuesta) %>%
  mutate(n_group = n()) %>%
  ungroup() %>%
  filter(n_group == 3) %>%
  select(-n_group)

# ------------------------------------------------------------------------------
# 6) Prepare factors to match Stata i. variables and ib2 reference
# + rescale age to reduce lme4 "nearly unidentifiable"
# ------------------------------------------------------------------------------
dt <- dt %>%
  mutate(
    id_f      = factor(idencuesta),
    clase_f   = relevel(factor(as.integer(clase)), ref = "2"),     # ib2
    ola_f     = factor(as.integer(ola), levels = c(1,3,6)),        # i.ola
    sexo_f    = factor(as.integer(m0_sexo)),                       # i.m0_sexo
    swb_f     = factor(as.integer(swb)),                           # i.swb
    emp_f     = factor(as.integer(employed)),                      # i.employed
    couple_f  = factor(as.integer(couple)),                        # i.couple
    tinst_f   = factor(as.integer(t02_01)),                        # i.t02_01

    # Correct label mapping: use clase (not class_num)
    class_lab = factor(as.integer(clase), levels = c(1,2,3),
                       labels = c("Clustering (\u03b2)",
                                  "Bridging (\u03b3)",
                                  "Isolated (\u03b1)")),

    # Rescaled age (centered decades)
    age10 = (as.numeric(m0_edad) - mean(as.numeric(m0_edad), na.rm = TRUE)) / 10,

    # For reverse (xtreg-like)
    clase_y = as.numeric(as.integer(clase))
  )

# Save balanced panel for debugging
write_csv(
  dt %>% select(idencuesta, ola, social_trust, trust_nb, clase, m0_sexo, m0_edad, age10, swb, employed, couple, t02_01),
  file.path(DIR_OUT, "balanced_panel_for_debug.csv")
)

# Diagnostics expected: N=3891
diag_lines <- c(
  "=== diagnostics (replicate_estimation_do_UPDATED.R) ===",
  paste0("Rows after merge: ", nrow(dt_merged)),
  paste0("Rows after filters+balance: ", nrow(dt)),
  paste0("Unique ids after balance: ", n_distinct(dt$idencuesta)),
  paste0("Waves present: ", paste(sort(unique(dt$ola)), collapse=", "))
)
writeLines(diag_lines, file.path(DIR_OUT, "diagnostics_estimation_do.txt"))
message(paste(diag_lines, collapse="\n"))

# ------------------------------------------------------------------------------
# 7) Fit models
#   Table 4: xtprobit (RE probit) + margins dydx(*)
#   Table A1: probit + margins dydx(*)
# ------------------------------------------------------------------------------
fit_re_probit <- function(outcome, controls = FALSE, nAGQ = 12L) {
  rhs <- if (!controls) {
    "clase_f + ola_f + (1 | id_f)"
  } else {
    "clase_f + sexo_f + age10 + swb_f + emp_f + couple_f + tinst_f + ola_f + (1 | id_f)"
  }
  fml <- as.formula(paste0(outcome, " ~ ", rhs))
  m <- glmer(
    fml, data = dt,
    family = binomial(link = "probit"),
    nAGQ = nAGQ,
    control = glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5))
  )
  vcv <- get_vcov_mer(m, cluster_var = "id_f")
  list(model = m, vcov = vcv, formula = deparse(fml))
}

fit_pool_probit <- function(outcome, controls = FALSE) {
  rhs <- if (!controls) {
    "clase_f + ola_f"
  } else {
    "clase_f + sexo_f + age10 + swb_f + emp_f + couple_f + tinst_f + ola_f"
  }
  fml <- as.formula(paste0(outcome, " ~ ", rhs))
  m <- glm(fml, data = dt, family = binomial(link = "probit"))
  vcv <- tryCatch(
    sandwich::vcovCL(m, cluster = dt$id_f, type = "HC1"),
    error = function(e) NULL
  )
  list(model = m, vcov = vcv, formula = deparse(fml))
}

# Table 4 (RE probit)
m1 <- fit_re_probit("social_trust", controls = FALSE, nAGQ = 12L)
m2 <- fit_re_probit("social_trust", controls = TRUE,  nAGQ = 12L)
m3 <- fit_re_probit("trust_nb",     controls = FALSE, nAGQ = 12L)
m4 <- fit_re_probit("trust_nb",     controls = TRUE,  nAGQ = 12L)

# Table A1 (pooled probit)
p1 <- fit_pool_probit("social_trust", controls = FALSE)
p2 <- fit_pool_probit("social_trust", controls = TRUE)
p3 <- fit_pool_probit("trust_nb",     controls = FALSE)
p4 <- fit_pool_probit("trust_nb",     controls = TRUE)

# ------------------------------------------------------------------------------
# 8) Marginal effects like Stata: margins, dydx(*)
# (Export only clase contrasts vs ref=2: "1 - 2" and "3 - 2")
# ------------------------------------------------------------------------------
extract_class_ames <- function(obj, outcome_name, model_name) {
  m <- obj$model
  v <- obj$vcov

  ame <- marginaleffects::avg_comparisons(
    m,
    variables = "clase_f",
    type = "response",
    vcov = v
  ) %>% as_tibble()

  ame %>%
    filter(grepl(" - 2$", contrast)) %>%
    mutate(
      outcome = outcome_name,
      model   = model_name,
      sig     = stars(p.value)
    ) %>%
    select(model, outcome, contrast, estimate, std.error, p.value, sig)
}

table4_R <- bind_rows(
  extract_class_ames(m1, "Generalized trust (social_trust)", "Model 1 (RE, no controls)"),
  extract_class_ames(m2, "Generalized trust (social_trust)", "Model 2 (RE, controls)"),
  extract_class_ames(m3, "Trust in neighbors (trust_nb)",    "Model 3 (RE, no controls)"),
  extract_class_ames(m4, "Trust in neighbors (trust_nb)",    "Model 4 (RE, controls)")
)
write_csv(table4_R, file.path(DIR_OUT, "table4_R.csv"))

tableA1_R <- bind_rows(
  extract_class_ames(p1, "Generalized trust (social_trust)", "Model 1 (Pooled, no controls)"),
  extract_class_ames(p2, "Generalized trust (social_trust)", "Model 2 (Pooled, controls)"),
  extract_class_ames(p3, "Trust in neighbors (trust_nb)",    "Model 3 (Pooled, no controls)"),
  extract_class_ames(p4, "Trust in neighbors (trust_nb)",    "Model 4 (Pooled, controls)")
)
write_csv(tableA1_R, file.path(DIR_OUT, "tableA1_R.csv"))

# ------------------------------------------------------------------------------
# 9) Figure 3 — predicted margins by class and wave (Stata: margins clase#ola; marginsplot, noci)
# Use the full-control RE models (m2, m4), like in Stata script
# ------------------------------------------------------------------------------
pred_grid <- function(obj, outcome_label) {
  m <- obj$model
  v <- obj$vcov
  pr <- marginaleffects::avg_predictions(
    m,
    by   = c("clase_f", "ola_f"),
    type = "response",
    vcov = v
  ) %>% as_tibble()

  pr %>%
    mutate(
      ola_raw = as.integer(as.character(ola_f)),
      wave_year = case_when(
        ola_raw == 1 ~ 2016L,
        ola_raw == 3 ~ 2018L,
        ola_raw == 6 ~ 2022L,
        TRUE ~ ola_raw
      ),
      class_num = as.integer(as.character(clase_f)),
      class_lab = factor(class_num, levels = c(1,2,3),
                         labels = c("Clustering (\u03b2)",
                                  "Bridging (\u03b3)",
                                  "Isolated (\u03b1)")),
      outcome_label = outcome_label
    )
}

preds <- bind_rows(
  pred_grid(m2, "Generalized trust (Pr=1)"),
  pred_grid(m4, "Trust in neighbors (Pr=1)")
)

fig3 <- ggplot(preds, aes(x = class_lab, y = estimate,
                          group = factor(wave_year), color = factor(wave_year))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  facet_wrap(~ outcome_label, scales = "free_y") +
  labs(
    x = "Class",
    y = "Predicted probability",
    color = NULL,
    title = ""
  ) +
  theme_ssr_big()



ggsave(file.path(DIR_OUT, "figure3_R.png"), fig3, width = 9, height = 5, dpi = 300)
ggsave(file.path(DIR_OUT, "figure3_R.pdf"), fig3, width = 9, height = 5)

# ------------------------------------------------------------------------------
# 10) TABLE A2 — Reverse relationship (trust -> class), like estimation.do
# Stata uses: xtreg clase i.social_trust ... i.ola, vce(cluster idencuesta)
# In R: lmer(clase_y ~ predictor + ... + ola_f + (1|id_f)) + CR1 cluster
# ------------------------------------------------------------------------------
fit_reverse_xtreg_like <- function(predictor, controls = FALSE) {
  # predictor is numeric 0/1 column name in dt; we treat it as factor (i.)
  pred_f <- paste0(predictor, "_f")
  d <- dt %>%
    mutate(
      social_trust_f = factor(social_trust, levels = c(0,1)),
      trust_nb_f     = factor(trust_nb,     levels = c(0,1))
    )

  rhs <- if (!controls) {
    paste0(pred_f, " + ola_f + (1 | id_f)")
  } else {
    paste0(pred_f, " + sexo_f + age10 + swb_f + emp_f + couple_f + tinst_f + ola_f + (1 | id_f)")
  }

  fml <- as.formula(paste0("clase_y ~ ", rhs))
  m <- lmer(
    fml, data = d,
    REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  vcv <- get_vcov_mer(m, cluster_var = "id_f")

  ame <- marginaleffects::avg_comparisons(
    m,
    variables = pred_f,
    vcov = vcv
  ) %>% as_tibble()

  list(model = m, vcov = vcv, ame = ame, formula = deparse(fml))
}

extract_a2 <- function(obj, model_name, predictor_label) {
  obj$ame %>%
    # keep the 1-0 comparison if present
    mutate(
      model = model_name,
      predictor = predictor_label,
      sig = stars(p.value)
    ) %>%
    select(model, predictor, term, contrast, estimate, std.error, p.value, sig)
}

ra1 <- fit_reverse_xtreg_like("social_trust", controls = FALSE)
ra2 <- fit_reverse_xtreg_like("social_trust", controls = TRUE)
ra3 <- fit_reverse_xtreg_like("trust_nb",     controls = FALSE)
ra4 <- fit_reverse_xtreg_like("trust_nb",     controls = TRUE)

tableA2_R <- bind_rows(
  extract_a2(ra1, "Model 1 (xtreg-like): clase ~ social_trust + i.ola", "Generalized trust (social_trust)"),
  extract_a2(ra2, "Model 2 (xtreg-like): + controls",                  "Generalized trust (social_trust)"),
  extract_a2(ra3, "Model 3 (xtreg-like): clase ~ trust_nb + i.ola",    "Trust in neighbors (trust_nb)"),
  extract_a2(ra4, "Model 4 (xtreg-like): + controls",                  "Trust in neighbors (trust_nb)")
)
write_csv(tableA2_R, file.path(DIR_OUT, "tableA2_R.csv"))

message("\nDONE: outputs written to /output")
message(" - table4_R.csv")
message(" - tableA1_R.csv")
message(" - tableA2_R.csv")
message(" - figure3_R.png/.pdf")
message(" - diagnostics_estimation_do.txt")