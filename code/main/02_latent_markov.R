# ── 02_latent_markov.R — Latent Markov model (LMest, K=3) ───────────────────
#
# Estimates K=1..5 LMMs; selects K=3 by BIC. Assigns alpha/beta/gamma labels
# by membership count (alpha = min) and domain entropy (gamma = max among rest).
# Anchors LMM sample to dt_analysis.rds to align with trust model sample.
#
# Outputs:
#   output/fit_table_cov.csv
#   output/latent_profiles_cov_K3.csv
#   output/transition_matrix_cov_K3.csv
#   output/classification_diagnostics_cov.txt
#   output/state_mapping_verification.txt
#   data/dt_states_cov.rds
#   data/posterior_probs_cov_std.rds

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(LMest)
  library(panelr)
  library(data.table)
  library(readr)
})

dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("data"),   showWarnings = FALSE, recursive = TRUE)

# Load canonical setup (STATE_MAP, helpers, etc.)
source(here::here("code", "00_setup.R"))

# ------------------------------------------------------------------------------
# USER KNOBS
# ------------------------------------------------------------------------------
K_GRID        <- 1:5
K_BASELINE    <- 3
SEED          <- 123
THRESH_STRICT <- 0.60    # strict assignment threshold on posterior max prob

# Suffix for output files: "" for main analysis, "_active_only" for sensitivity S7
OUT_SUFFIX <- if (exists("MEMBER_CODE_LOGIC") && MEMBER_CODE_LOGIC == "active_only") "_active_only" else ""
message("  02_latent_markov.R | MEMBER_CODE_LOGIC = '", MEMBER_CODE_LOGIC,
        "' | output suffix = '", OUT_SUFFIX, "'")

# MAP_GAMMA_BY: method for identifying the bridging state.
#   "entropy"         = gamma is the state with highest domain entropy (most diverse)
#   "mean_membership" = gamma is the state with the highest mean membership count
# Recommended: "entropy" captures structural diversity, not just quantity.
MAP_GAMMA_BY <- "entropy"

# ── 0. Load data ────────────────────────────────────────────────────────────
stopifnot(file.exists(here::here("data/ELSOC_Long.RData")))
obj_names <- load(here::here("data/ELSOC_Long.RData"))

a <- if ("elsoc_long_2016_2022" %in% obj_names) {
  get("elsoc_long_2016_2022")
} else if (length(obj_names) == 1) {
  get(obj_names[1])
} else {
  stop("Cannot identify ELSOC object. Available: ",
       paste(obj_names, collapse = ", "))
}

# ── 0b. Sample filters (must match 01_descriptive_stats.R) ─────────────────
# Without these filters the LMM sample diverges from the probit sample,
# producing ~270 extra individuals without valid trust data.
n_raw <- nrow(a)
if ("muestra" %in% names(a)) {
  a <- a %>% filter(muestra == 1)
  message("Filter muestra==1: ", n_raw, " -> ", nrow(a), " obs")
} else {
  warning("Variable 'muestra' not found in ELSOC_Long. Filter not applied.")
}
n_after_muestra <- nrow(a)
if ("tipo_atricion" %in% names(a)) {
  a <- a %>% filter(tipo_atricion == 1)
  message("Filter tipo_atricion==1: ", n_after_muestra, " -> ", nrow(a), " obs")
} else {
  warning("Variable 'tipo_atricion' not found in ELSOC_Long. Filter not applied.")
}

# ── 1. Prepare panel (waves 1,3,6 -> 1,2,3) ─────────────────────────────────
# Missing codes recoded to NA before binary c12 recode to prevent collapse.
a_sub <- a %>%
  filter(ola %in% c(1, 3, 6)) %>%
  mutate(
    ola = case_when(ola == 1 ~ 1L,
                    ola == 3 ~ 2L,
                    ola == 6 ~ 3L)
  ) %>%
  select(
    idencuesta, ola, ponderador02,
    m0_sexo, m0_edad, m01,
    c02, c03, c04,
    c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08
  )

# Missing codes -> NA before binary recode
a_sub <- a_sub %>%
  mutate(across(where(is.numeric),
                ~ replace(.x, .x %in% c(-999, -888, -777, -666), NA)))

a_full <- a_sub %>%
  # Use member_binary() from 00_setup.R to respect MEMBER_CODE_LOGIC
  # (any_member: c12 >= 2; active_only: c12 == 3 — see sensitivity S7)
  mutate(across(matches("^c12"), ~ member_binary(.x))) %>%
  mutate(conf_gral = case_when(c02 == 1 ~ 1L,
                               c02 == 2 ~ 0L,
                               c02 == 3 ~ 1L,
                               TRUE ~ NA_integer_)) %>%
  mutate(
    mujer = case_when(m0_sexo == 1 ~ 0L,
                      m0_sexo == 2 ~ 1L,
                      TRUE ~ NA_integer_),
    edad = case_when(m0_edad %in% 18:24  ~ "18_24",
                     m0_edad %in% 25:34  ~ "25_34",
                     m0_edad %in% 35:44  ~ "35_44",
                     m0_edad %in% 45:54  ~ "45_54",
                     m0_edad %in% 55:64  ~ "55_64",
                     m0_edad %in% 65:120 ~ "65",
                     TRUE ~ NA_character_),
    nivel_educ = case_when(m01 %in% 1:3   ~ "básica",
                           m01 %in% 4:5   ~ "media",
                           m01 %in% 6:7   ~ "técnica",
                           m01 %in% 8:10  ~ "univers",
                           TRUE ~ NA_character_)
  ) %>%
  select(
    idencuesta, ola, ponderador02,
    mujer, edad, nivel_educ,
    c02, c03, c04, conf_gral,
    c12_01, c12_02, c12_03, c12_04, c12_05, c12_06, c12_07, c12_08
  )

# ── 1b. Anchor LMM sample to dt_analysis.rds ─────────────────────────────────
# Restricts a_full to IDs that passed all filters in 01_descriptive_stats.R,
# ensuring the LMM and trust models operate on the same N=1,297 individuals.
dt_analysis_path <- here::here("data", "dt_analysis.rds")
if (file.exists(dt_analysis_path)) {
  dt_anal_ids <- readRDS(dt_analysis_path) %>%
    dplyr::distinct(id) %>%
    dplyr::mutate(idencuesta = as.character(id)) %>%
    dplyr::pull(idencuesta)

  n_before_anchor <- length(unique(a_full$idencuesta))
  a_full <- a_full %>%
    dplyr::mutate(idencuesta_chr = as.character(idencuesta)) %>%
    dplyr::filter(idencuesta_chr %in% dt_anal_ids) %>%
    dplyr::select(-idencuesta_chr)
  n_after_anchor <- length(unique(a_full$idencuesta))

  message("Anchoring to dt_analysis.rds: ",
          n_before_anchor, " -> ", n_after_anchor, " individuos")
  message("  (", n_before_anchor - n_after_anchor,
          " individuals without valid trust/controls excluded from LMM)")
} else {
  warning("dt_analysis.rds not found. ",
          "Run 01_descriptive_stats.R first to anchor the LMM sample.")
  # Fallback: drop_na on LMM covariates (may produce incorrect N)
  a_full <- a_full %>% drop_na(mujer, edad, nivel_educ)
}

# ── 2. Enforce balanced panel (3 waves) ───────────────────────────────────────
items_check <- c("c12_01","c12_02","c12_03","c12_04",
                  "c12_05","c12_06","c12_07","c12_08")
a_full <- panel_data(a_full, id = idencuesta, wave = ola) %>%
  complete_data(min.waves = "all", vars = items_check) %>%
  as.data.frame()

# lmestSearch() requires no NAs in transition covariates.
# Drop observations with NA in mujer/edad/nivel_educ and re-enforce balance.
n_pre_na <- length(unique(a_full$idencuesta))
a_full <- a_full %>%
  dplyr::filter(!is.na(mujer), !is.na(edad), !is.na(nivel_educ))
# Re-enforce balance: keep only IDs with exactly 3 waves after the drop
ids_balanced <- a_full %>%
  dplyr::count(idencuesta) %>%
  dplyr::filter(n == 3) %>%
  dplyr::pull(idencuesta)
a_full <- a_full %>% dplyr::filter(idencuesta %in% ids_balanced)
n_post_na <- length(unique(a_full$idencuesta))
if (n_pre_na - n_post_na > 0)
  message("  drop_na(LMM covariates): ", n_pre_na, " -> ", n_post_na,
          " ids (", n_pre_na - n_post_na, " excluded due to NA in mujer/edad/nivel_educ)")

message("Balanced panel: ", nrow(a_full), " rows | ",
        length(unique(a_full$idencuesta)), " ids | ",
        length(unique(a_full$ola)), " waves")
message("  Target (Stata probit): ~3,891 obs (1,297 individuals x 3 waves)")
message("  If N differs substantially, check that 01 used the same filters as estimation.do")

# ── 3. Latent Markov model selection (K=1..5) ─────────────────────────────────
items <- c("c12_01","c12_02","c12_03","c12_04","c12_05","c12_06","c12_07","c12_08")

stopifnot(all(items %in% names(a_full)))
stopifnot(all(c("mujer","edad","nivel_educ","idencuesta","ola") %in% names(a_full)))

set.seed(SEED)
mod_sel <- lmestSearch(
  responsesFormula = c12_01 + c12_02 + c12_03 + c12_04 + c12_05 + c12_06 + c12_07 + c12_08 ~ NULL,
  latentFormula    = ~ mujer + edad + nivel_educ,
  index            = c("idencuesta", "ola"),
  data             = a_full,
  output           = TRUE,
  out_se           = TRUE,
  version          = "categorical",
  paramLatent      = "multilogit",
  k                = K_GRID,
  seed             = SEED
)

fit_stats <- tibble(
  K      = K_GRID,
  logLik = as.numeric(mod_sel[["lkv"]]),
  AIC    = as.numeric(mod_sel[["Aic"]]),
  BIC    = as.numeric(mod_sel[["Bic"]])
)

np_try <- mod_sel[["np"]]
fit_stats <- fit_stats %>%
  mutate(npar = if (!is.null(np_try)) as.numeric(np_try) else NA_real_)

write_csv(fit_stats, here::here("output", paste0("fit_table_cov", OUT_SUFFIX, ".csv")))
print(fit_stats)

K_star <- fit_stats$K[which.min(fit_stats$BIC)]
message("Selected K by BIC: K* = ", K_star, " | Baseline K=", K_BASELINE)

# ── 4. Extract K=3 model ──────────────────────────────────────────────────────
if (!is.null(mod_sel$out.single) && length(mod_sel$out.single) >= K_BASELINE) {
  modeloK <- mod_sel$out.single[[K_BASELINE]]
} else {
  stop("mod_sel$out.single not found or does not contain K=", K_BASELINE)
}

Psi <- modeloK$Psi
Pi  <- modeloK$Pi
V   <- modeloK$V

# ── 5. Standardize posterior array V to n x K x TT ───────────────────────────
reshape_V <- function(V, n, k, TT) {
  d <- dim(V)
  if (length(d) != 3) stop("dim(V) not 3D: ", paste(d, collapse=" x "))

  perms <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  best <- NULL
  best_err <- Inf
  for (p in perms) {
    Vp <- aperm(V, p)
    if (!all(dim(Vp) == c(n,k,TT))) next
    err <- max(abs(apply(Vp, c(1,3), sum) - 1), na.rm = TRUE)
    if (err < best_err) { best_err <- err; best <- Vp }
  }
  if (is.null(best)) stop("Could not reshape V to (n,k,TT). dims=", paste(d, collapse="x"))
  best
}

ids   <- sort(unique(a_full$idencuesta))
waves <- sort(unique(a_full$ola))
n     <- length(ids)
TT    <- length(waves)
k     <- K_BASELINE
r     <- length(items)

V_std <- reshape_V(V, n=n, k=k, TT=TT)

# ── 6. Extract item-response probabilities P(Y=1|state) ───────────────────────
perm_try <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
Psi_std <- NULL
for (p in perm_try) {
  Pp <- aperm(Psi, p)
  if (all(dim(Pp) == c(r,k,2))) { Psi_std <- Pp; break }
}
if (is.null(Psi_std)) stop("Could not reshape Psi to (r,k,2). dims=",
                            paste(dim(Psi), collapse="x"))

prob1 <- Psi_std[,,2, drop=TRUE]  # r x k, probability of category "1"
rownames(prob1) <- items

# ── 7. Map states to alpha/beta/gamma ─────────────────────────────────────────
# alpha = min expected membership count (colSums lowest)
# gamma = max domain entropy among remaining states (most diverse)
# beta  = residual state

# Expected membership count per state
exp_count <- colSums(prob1, na.rm = TRUE)
names(exp_count) <- paste0("State_", 1:k)

# Normalised domain entropy per state
state_entropy <- sapply(1:k, function(s) {
  p <- prob1[, s]
  p <- p / (sum(p) + 1e-12)
  -sum(p * log(p + 1e-12)) / log(length(p))
})
names(state_entropy) <- paste0("State_", 1:k)

# alpha = lowest count
alpha_state <- which.min(exp_count)

# gamma = highest entropy among non-alpha states
remaining <- setdiff(1:k, alpha_state)
gamma_state <- remaining[which.max(state_entropy[remaining])]

# beta = residual state
beta_state <- setdiff(1:k, c(alpha_state, gamma_state))

# Guard: beta should be scalar
if (length(beta_state) != 1) {
  warning("beta_state tiene longitud inesperada: ", length(beta_state),
          ". Usando order(exp_count) como fallback.")
  ord         <- order(exp_count)
  alpha_state <- ord[1]
  gamma_state <- ord[3]
  beta_state  <- ord[2]
}

message("\n=== STATE MAPPING ===")
message("  \u03b1 (Isolate)  = State ", alpha_state,
        " | exp_count=", round(exp_count[alpha_state], 3),
        " | entropy=",   round(state_entropy[alpha_state], 3))
message("  \u03b2 (Closed)   = State ", beta_state,
        " | exp_count=", round(exp_count[beta_state], 3),
        " | entropy=",   round(state_entropy[beta_state], 3))
message("  \u03b3 (Bridging) = State ", gamma_state,
        " | exp_count=", round(exp_count[gamma_state], 3),
        " | entropy=",   round(state_entropy[gamma_state], 3))

# Verify that gamma has high prob1 in >= 4 domains
gamma_n_high <- sum(prob1[, gamma_state] >= 0.30)
beta_n_high  <- sum(prob1[, beta_state]  >= 0.30)
alpha_n_high <- sum(prob1[, alpha_state] >= 0.30)

message("\n  Verification (prob1 >= 0.30 per state):")
message("    \u03b1: ", alpha_n_high, " active domains (expected: <= 2)")
message("    \u03b2: ", beta_n_high,  " active domains (expected: 1-3, concentrated)")
message("    \u03b3: ", gamma_n_high, " active domains (expected: >= 4, diverse)")

if (gamma_n_high < 3) {
  warning("gamma_state has only ", gamma_n_high,
          " domains with prob1>=0.30. Check that the mapping is correct.",
          "\nConsider changing MAP_GAMMA_BY or K_BASELINE.")
}
if (alpha_n_high > 3) {
  warning("alpha_state has ", alpha_n_high,
          " active domains -- may not correspond to Isolate.")
}

# ── 8. Build dt_states_cov with posteriors and modal assignments ────────────
Uhat <- sapply(seq_len(TT), function(t)
  max.col(V_std[, , t, drop = TRUE], ties.method = "first"))
Uhat <- matrix(Uhat, nrow = n, ncol = TT)

label_state <- function(s) {
  dplyr::case_when(
    s == alpha_state ~ "alpha",
    s == gamma_state ~ "gamma",
    TRUE             ~ "beta"
  )
}

# Soft posteriors (long format)
post_long <- rbindlist(lapply(seq_len(TT), function(t) {
  pp <- V_std[, , t, drop = TRUE]
  data.table(
    idencuesta = ids,
    ola        = waves[t],
    p_alpha    = pp[, alpha_state],
    p_beta     = pp[, beta_state],
    p_gamma    = pp[, gamma_state],
    p_max      = apply(pp, 1, max)
  )
}))

# Modal (hard) assignments (long format)
modal_dt <- rbindlist(lapply(seq_len(TT), function(t) {
  data.table(
    idencuesta = ids,
    ola        = waves[t],
    state_num  = Uhat[, t],
    position   = label_state(Uhat[, t])
  )
}))

dt_states <- as.data.table(a_full)
dt_states <- merge(dt_states, post_long, by = c("idencuesta","ola"), all.x = TRUE)
dt_states <- merge(dt_states, modal_dt,  by = c("idencuesta","ola"), all.x = TRUE)
dt_states[, position := factor(position, levels = c("alpha","beta","gamma"))]

# Strict classification: NA if max posterior < THRESH_STRICT
dt_states[, position_strict := as.character(position)]
dt_states[p_max < THRESH_STRICT, position_strict := NA_character_]
dt_states[, position_strict := factor(position_strict, levels = c("alpha","beta","gamma"))]

# ── 9. Diagnostics and exports ────────────────────────────────────────────────
pp_mat <- as.matrix(dt_states[, .(p_alpha, p_beta, p_gamma)])
eps <- 1e-12
entropy_norm <- -mean(rowSums(pp_mat * log(pp_mat + eps), na.rm = TRUE)) / log(k)

avg_pp_tbl <- rbindlist(lapply(levels(dt_states$position), function(pos) {
  col_nm <- paste0("p_", pos)
  sub <- dt_states[position == pos]
  data.table(
    position   = pos,
    n_assigned = nrow(sub),
    avg_pp     = mean(sub[[col_nm]], na.rm = TRUE),
    pct_total  = nrow(sub) / nrow(dt_states) * 100
  )
}))

quality_global <- dt_states[, .(
  mean_p_max   = mean(p_max, na.rm = TRUE),
  share_strict = mean(!is.na(position_strict)),
  n            = .N
)]

quality_by_wave <- dt_states[, .(
  mean_p_max   = mean(p_max, na.rm = TRUE),
  share_strict = mean(!is.na(position_strict)),
  n            = .N
), by = ola][order(ola)]

# Average transition matrix
Pi_avg <- Pi
if (length(dim(Pi)) == 3) Pi_avg <- apply(Pi, c(1,2), mean)

ord_states <- c(alpha_state, beta_state, gamma_state)
Pi_reord <- Pi_avg[ord_states, ord_states, drop = FALSE]
rownames(Pi_reord) <- colnames(Pi_reord) <- c("alpha","beta","gamma")

# Save item-response profiles and transition matrix
prob_df <- as.data.frame(prob1)
colnames(prob_df) <- paste0("State_", 1:k)
prob_df$item       <- items
write_csv(prob_df, here::here("output", paste0("latent_profiles_cov_K3", OUT_SUFFIX, ".csv")))

trans_df <- as.data.frame(Pi_reord) %>% tibble::rownames_to_column("from")
write_csv(trans_df, here::here("output", paste0("transition_matrix_cov_K3", OUT_SUFFIX, ".csv")))

# ── 9b. State mapping verification report ─────────────────────────────────────
# Review output/state_mapping_verification.txt before running 03_trust_models.R

prob_display <- round(prob1, 3)
rownames(prob_display) <- items
colnames(prob_display) <- c(paste0("State_", 1:k))

verify_lines <- c(
  "=== STATE MAPPING VERIFICATION (002_long_latent_class.R) ===",
  paste0("Dataset: ELSOC Long 2016-2022 | K=", K_BASELINE, " | MAP_GAMMA_BY=", MAP_GAMMA_BY),
  "",
  "--- Mapping (STATE_MAP canonical) ---",
  paste0("  \u03b1 (Isolate/Apathetic)  -> LMM State ", alpha_state, "  | clase=3 | ref=ib3"),
  paste0("  \u03b2 (Closed/Clustering)  -> LMM State ", beta_state,  "  | clase=1 | ref=ib1"),
  paste0("  \u03b3 (Bridging/Broker)    -> LMM State ", gamma_state, "  | clase=2 | ref=ib2 [REFERENCIA]"),
  "",
  "--- Descriptivos por estado ---",
  paste0("  exp_count  (colSums prob1): ", paste(round(exp_count, 3), collapse=" | ")),
  paste0("  entropy    (domain divers): ", paste(round(state_entropy, 3), collapse=" | ")),
  "",
  "--- Membership probabilities P(Y=1|state) ---",
  "    Items x States:",
  capture.output(print(prob_display)),
  "",
  "--- Dominios activos (prob1 >= 0.30) ---",
  paste0("  \u03b1 State_", alpha_state, ": ", alpha_n_high, " active domains (expected: <= 2)"),
  paste0("  \u03b2 State_", beta_state,  ": ", beta_n_high,  " active domains (expected: 1-3, concentrated)"),
  paste0("  \u03b3 State_", gamma_state, ": ", gamma_n_high, " active domains (expected: >= 4, diverse)"),
  "",
  "--- Substantive check (answer Yes/No before running 03_trust_models.R) ---",
  "  1. Does alpha have low probabilities across most domains?",
  "  2. Does beta concentrate memberships in 1-3 domains (e.g., religious, neighbourhood)?",
  "  3. Does gamma have probabilities >= 0.30 in >= 4 distinct domains?",
  "  4. Does the alpha/beta/gamma ordering make sense as Isolate/Closed/Bridging?",
  "",
  "  If NO: change MAP_GAMMA_BY, review K_BASELINE, or adjust alpha_state manually.",
  "",
  "--- Average transition matrix (alpha/beta/gamma) ---",
  capture.output(print(round(Pi_reord, 3))),
  "",
  "--- Classification quality ---",
  capture.output(print(quality_global)),
  paste0("  Entropy norm (lower=better): ", round(entropy_norm, 4)),
  paste0("  Strict threshold (p_max >= ", THRESH_STRICT, "): ",
         round(quality_global$share_strict, 3), " obs clasificados")
)
writeLines(verify_lines, here::here("output", paste0("state_mapping_verification", OUT_SUFFIX, ".txt")))
message("\n  Verification saved: output/state_mapping_verification.txt")
message("  *** Review this file before running 03_trust_models.R ***")

# Diagnostics completo
diag_txt <- c(
  "=== Covariate LMM diagnostics (K=3) ===",
  paste0("Balanced panel: ids=", n, " | waves=", TT, " | items=", r),
  paste0("Normalized entropy (lower=better) = ", round(entropy_norm, 4)),
  paste0("Strict threshold p_max >= ", THRESH_STRICT),
  paste0("MAP_GAMMA_BY = ", MAP_GAMMA_BY),
  "",
  "State mapping (canonical STATE_MAP):",
  paste0("  \u03b1 (Isolate)  = State ", alpha_state),
  paste0("  \u03b2 (Closed)   = State ", beta_state),
  paste0("  \u03b3 (Bridging) = State ", gamma_state, " [REFERENCIA ib2]"),
  "",
  "Assignment quality (by assigned class):",
  capture.output(print(avg_pp_tbl)),
  "",
  "Quality (global):",
  capture.output(print(quality_global)),
  "",
  "Quality by wave:",
  capture.output(print(quality_by_wave)),
  "",
  "Average transition matrix (alpha/beta/gamma):",
  capture.output(print(round(Pi_reord, 3)))
)
writeLines(diag_txt, here::here("output", paste0("classification_diagnostics_cov", OUT_SUFFIX, ".txt")))

# ── 10. Save objects ──────────────────────────────────────────────────────────
saveRDS(as_tibble(dt_states), here::here("data", paste0("dt_states_cov", OUT_SUFFIX, ".rds")))

saveRDS(
  list(
    fit_table      = fit_stats,
    mod_sel        = mod_sel,
    modelo3        = modeloK,
    ids            = ids,
    waves          = waves,
    V_std          = V_std,
    Psi            = Psi,
    Pi             = Pi,
    prob1          = prob1,
    prob_df        = prob_df,
    Pi_reord       = Pi_reord,
    exp_count      = exp_count,
    state_entropy  = state_entropy,
    alpha_state    = alpha_state,
    beta_state     = beta_state,
    gamma_state    = gamma_state,
    alpha_n_high   = alpha_n_high,
    beta_n_high    = beta_n_high,
    gamma_n_high   = gamma_n_high,
    entropy_norm   = entropy_norm,
    quality_global = quality_global,
    quality_by_wave = quality_by_wave,
    THRESH_STRICT  = THRESH_STRICT,
    MAP_GAMMA_BY   = MAP_GAMMA_BY,
    MEMBER_CODE_LOGIC = MEMBER_CODE_LOGIC,
    STATE_MAP      = STATE_MAP
  ),
  here::here("data", paste0("posterior_probs_cov_std", OUT_SUFFIX, ".rds"))
)

message("\n[02_latent_markov.R] done. MEMBER_CODE_LOGIC = '", MEMBER_CODE_LOGIC, "'")
message("Analysis N target: 1,297 individuals, 3,891 person-waves")
message("Saved: output/fit_table_cov",               OUT_SUFFIX, ".csv")
message("Saved: output/latent_profiles_cov_K3",      OUT_SUFFIX, ".csv")
message("Saved: output/transition_matrix_cov_K3",    OUT_SUFFIX, ".csv")
message("Saved: output/classification_diagnostics_cov", OUT_SUFFIX, ".txt")
message("Saved: output/state_mapping_verification",  OUT_SUFFIX, ".txt")
message("Saved: data/dt_states_cov",                 OUT_SUFFIX, ".rds")
message("Saved: data/posterior_probs_cov_std",       OUT_SUFFIX, ".rds")
message("")
message("Entropy=",   round(entropy_norm, 4),
        " | mean_p_max=",    round(quality_global$mean_p_max, 3),
        " | share_strict=",  round(quality_global$share_strict, 3))
message("")
message("Final mapping (confirm in state_mapping_verification.txt):")
message("  \u03b1 (Isolate)  = State ", alpha_state, " | clase=3")
message("  \u03b2 (Closed)   = State ", beta_state,  " | clase=1")
message("  \u03b3 (Bridging) = State ", gamma_state, " | clase=2 [REFERENCIA]")




# ── 11. Figure 1 — Latent class profiles ─────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
  library(glue)
})

# Domain labels for figure
domain_order_fig <- c("c12_01","c12_02","c12_03","c12_04",
                       "c12_05","c12_06","c12_07","c12_08")

domain_labels_fig <- c(
  c12_01 = "Neighborhood orgs.",
  c12_02 = "Religious orgs.",
  c12_03 = "Political parties",
  c12_04 = "Labor unions",
  c12_05 = "Professional assoc.",
  c12_06 = "Charitable orgs.",
  c12_07 = "Sports clubs",
  c12_08 = "Student orgs."
)

# State key: maps section-7 indices to labels and colors
state_key_fig <- tibble(
  state_idx  = c(alpha_state,   beta_state,      gamma_state),
  state_col  = paste0("State_", c(alpha_state, beta_state, gamma_state)),
  state_long = c("\u03b1 (isolation)",
                 "\u03b2 (clustering)",
                 "\u03b3 (bridging)"),
  col        = c("#4E79A7", "#59A14F", "#E15759")
)

state_pal_fig  <- setNames(state_key_fig$col, state_key_fig$state_long)
state_lvls_fig <- state_key_fig$state_long

# Prevalence from dt_states (used for strip subtitles)
prev_tbl <- tryCatch({
  readRDS(here::here("data","dt_states_cov.rds")) %>%
    count(position) %>%
    mutate(
      pct = n / sum(n),
      state_long = case_when(
        position == "alpha" ~ "\u03b1 (isolation)",
        position == "beta"  ~ "\u03b2 (clustering)",
        position == "gamma" ~ "\u03b3 (bridging)"
      )
    )
}, error = function(e) NULL)

make_strip_label <- function(sl) {
  if (!is.null(prev_tbl)) {
    row <- prev_tbl %>% filter(state_long == sl)
    if (nrow(row) == 1)
      return(glue("{sl}  [{scales::percent(row$pct, accuracy = 1)}]"))
  }
  sl
}
state_labels_fig <- setNames(
  sapply(state_lvls_fig, make_strip_label),
  state_lvls_fig
)

# Reshape prob1 to long for plotting
prob_long_fig <- as.data.frame(prob1) %>%
  setNames(paste0("State_", seq_len(ncol(prob1)))) %>%
  mutate(item = rownames(prob1)) %>%
  pivot_longer(starts_with("State_"),
               names_to  = "state_col",
               values_to = "prob") %>%
  inner_join(state_key_fig, by = "state_col") %>%
  mutate(
    state_long = factor(state_long, levels = state_lvls_fig),
    item_lab   = factor(
      recode(item, !!!domain_labels_fig),
      levels = rev(recode(domain_order_fig, !!!domain_labels_fig))
    ),
    prob_label = if_else(prob >= 0.15,
                         percent(prob, accuracy = 1),
                         NA_character_)
  )

# Plot
p_profiles_fig <- ggplot(prob_long_fig,
                         aes(x = item_lab, y = prob,
                             color = state_long, fill = state_long)) +

  # Subtle bar fill behind lollipop (same color, very transparent)
  #geom_col(width = 0.55, alpha = 0.08, show.legend = FALSE) +

  # Lollipop stem
  geom_segment(aes(xend = item_lab, y = 0, yend = prob),
             linewidth = 0.5,      # ← antes era 1.1
             alpha = 0.70,
             lineend = "round", show.legend = FALSE) +

  # Reference dashed line at 50%
  geom_hline(yintercept = 0.50, linetype = "22",
             color = "grey55", linewidth = 0.5) +

  # Dot
  geom_point(size = 4.5, alpha = 0.95, show.legend = FALSE) +

  # Percentage label above dot (only for prob >= 15%)
  geom_text(aes(label = prob_label, y = prob + 0.065),
            size = 2.9, fontface = "bold",
            na.rm = TRUE, show.legend = FALSE) +

  coord_flip() +

  facet_wrap(~ state_long, ncol = 1, scales = "fixed",
             labeller = labeller(state_long = state_labels_fig)) +

  scale_color_manual(values = state_pal_fig) +
  scale_fill_manual( values = state_pal_fig) +

  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, .25, .50, .75, 1.0),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.10))
  ) +

  labs(
    x = NULL,
    y = "Pr(member\u202f|\u202fstate)",
    caption = paste0(
      "Item-response probabilities from Latent Markov Model (K\u2009=\u2009",
      K_BASELINE, "; ELSOC balanced panel 2016\u20132022). ",
      "Labels shown for Pr\u202f\u2265\u202f15%. Dashed line at 50%.\n",
      "Strip brackets show pooled position prevalence. ",
      "\u03b1\u202f=\u202fisolation, \u03b2\u202f=\u202fclustering, \u03b3\u202f=\u202fbridging."
    )
  ) +

  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    strip.text         = element_text(face = "bold", size = 12.5,
                                      color = "white",
                                      margin = margin(t = 7, b = 7)),
    strip.background   = element_rect(color = NA, fill = "grey30"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    panel.spacing      = unit(1.2, "lines"),
    axis.text.y        = element_text(size = 11, color = "grey15"),
    axis.text.x        = element_text(size = 10, color = "grey30"),
    axis.title.x       = element_text(face = "bold", size = 12,
                                      margin = margin(t = 10)),
    plot.caption       = element_text(size = 8, color = "grey50",
                                      hjust = 0, lineheight = 1.3,
                                      margin = margin(t = 12)),
    plot.margin        = margin(t = 10, r = 22, b = 10, l = 10)
  )

# Colour strip backgrounds to match state palette
gt_profiles <- tryCatch({
  gb <- ggplot_build(p_profiles_fig)
  gt <- ggplot_gtable(gb)
  strip_idx <- which(grepl("^strip", gt$layout$name))

  for (si in strip_idx) {
    # Extract text label from strip grob tree
    strip_label <- tryCatch(
      gt$grobs[[si]]$grobs[[1]]$children[[2]]$children[[1]]$label,
      error = function(e) ""
    )
    # Match by checking which state_lvl is a prefix of the strip label
    matched <- state_lvls_fig[
      sapply(state_lvls_fig, function(s) startsWith(strip_label, s))
    ]
    if (length(matched) == 1) {
      gt$grobs[[si]]$grobs[[1]]$children[[1]]$gp$fill <- state_pal_fig[[matched]]
    }
  }
  gt
}, error = function(e) {
  message("  Strip color override failed: ", conditionMessage(e))
  ggplot_gtable(ggplot_build(p_profiles_fig))
})

# Save
png(here::here("output", "fig_profiles_clean.png"),
    width = 7.0, height = 10.0, units = "in", res = 300)
grid::grid.draw(gt_profiles)
dev.off()

pdf(here::here("output", "fig_profiles_clean.pdf"),
    width = 7.0, height = 10.0)
grid::grid.draw(gt_profiles)
dev.off()

message("  Figure 1 saved: output/fig_profiles_clean.png/.pdf")
