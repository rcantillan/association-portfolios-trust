# ==============================================================================
# 02_latent_markov.R — Latent Markov model for membership portfolios (α/β/γ)
# ------------------------------------------------------------------------------
# FIXES INCLUDED (to make everything run and avoid "all alpha" bug):
#   (1) Loads all needed packages explicitly (incl. data.table) so .() works.
#   (2) Robust reconstruction of modal states from posterior probs V when k == TT.
#       Uses "sums to 1 across states" diagnostic to identify state dimension.
#   (3) Uses max.col() on a proper 2D matrix (drop=TRUE), not a 3D array.
#   (4) Guardrails: checks that position is not degenerate; prints counts by position.
# ------------------------------------------------------------------------------

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(LMest)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(readr)
  library(data.table)
})

# ------------------------------------------------------------------------------
# Parameters (edit if needed)
# ------------------------------------------------------------------------------
K_GRID      <- 1:5
K_BASELINE  <- 3
N_STARTS    <- 30     # multi-start per K (recommended for mixture models)
MOD_TRANS   <- 1      # 1 = homogeneous transitions over time; 0 = time-varying
TOL         <- 1e-8
MAXIT       <- 1000

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
load_elsoc_object <- function(path) {
  stop_if_missing(c(path))
  obj_names <- load(path)
  if ("elsoc_long_2016_2022" %in% obj_names) return(get("elsoc_long_2016_2022"))
  if (length(obj_names) == 1) return(get(obj_names[1]))
  stop("Could not find elsoc_long_2016_2022 in ", path, ". Objects: ", paste(obj_names, collapse = ", "))
}

fit_best_lmest <- function(S, yv, k, nstart = 20, mod = 1, tol = 1e-8, maxit = 1000) {
  best <- NULL
  best_lk <- -Inf
  for (s in seq_len(nstart)) {
    fit <- LMest::est_lm_basic(S, yv, k, start = 1, mod = mod, tol = tol, maxit = maxit, out_se = FALSE)
    if (is.finite(fit$lk) && fit$lk > best_lk) {
      best <- fit
      best_lk <- fit$lk
    }
  }
  if (is.null(best)) stop("LMest failed for k=", k, " (no finite logLik across starts).")
  best
}

# Robustly reshape V to n x k x TT using "sum to 1 across states" diagnostic
reshape_V_to_nkT <- function(V, n, k, TT) {
  dV <- dim(V)
  if (length(dV) != 3) stop("Unexpected dim(V) = ", paste(dV, collapse=" x "), ". Expected 3D array.")

  # Put respondent dimension first
  dim_n <- if (any(dV == n)) which(dV == n)[1] else which.max(dV)
  Vn <- aperm(V, c(dim_n, setdiff(1:3, dim_n)))  # n x d2 x d3

  d2 <- dim(Vn)[2]
  d3 <- dim(Vn)[3]

  # Candidate A: states are dim2
  sumA <- apply(Vn, c(1,3), sum)                 # n x d3
  errA <- max(abs(sumA - 1), na.rm = TRUE)

  # Candidate B: states are dim3
  sumB <- apply(Vn, c(1,2), sum)                 # n x d2
  errB <- max(abs(sumB - 1), na.rm = TRUE)

  if (errA <= errB) {
    # Vn is n x (states) x (time) OR n x (time) x (states) but states identified as dim2
    V_std <- Vn
    # If lengths look swapped (rare but possible), align to n x k x TT
    if (!(dim(V_std)[2] == k && dim(V_std)[3] == TT) && (dim(V_std)[2] == TT && dim(V_std)[3] == k)) {
      V_std <- aperm(Vn, c(1,3,2))
    }
  } else {
    # states are dim3 => Vn is n x (time) x (states) => permute to n x k x TT
    V_std <- aperm(Vn, c(1,3,2))
  }

  # Guardrail: sums to 1 across states (dim2) for each (i,t)
  check <- apply(V_std, c(1,3), sum)
  if (max(abs(check - 1), na.rm = TRUE) > 1e-3) {
    stop("Posterior probs do not sum to 1 across states after reshaping. Check V dimensions.")
  }

  # Additional size sanity (not fatal if TT/k differ in edge cases, but should match here)
  if (!(dim(V_std)[1] == n && dim(V_std)[2] == k && dim(V_std)[3] == TT)) {
    warning("V_std dims are ", paste(dim(V_std), collapse=" x "),
            " but expected ", paste(c(n,k,TT), collapse=" x "),
            ". Proceeding because sum-to-1 check passed.")
  }

  V_std
}

# ------------------------------------------------------------------------------
# 1) Load raw data
# ------------------------------------------------------------------------------
raw <- load_elsoc_object(here::here("data", "ELSOC_Long.RData"))
dt_raw <- dplyr::as_tibble(raw)

# ------------------------------------------------------------------------------
# 2) Restrict to analysis waves (ola 1,3,6 -> 1,2,3) and standardize id
# ------------------------------------------------------------------------------
dt <- dt_raw %>%
  dplyr::filter(ola %in% WAVES_RAW) %>%
  dplyr::mutate(
    ola = dplyr::case_when(ola == 1 ~ 1L, ola == 3 ~ 2L, ola == 6 ~ 3L),
    id  = idencuesta
  )

# ------------------------------------------------------------------------------
# 3) Membership domains (8) from c12_01..c12_08
# ------------------------------------------------------------------------------
missing_member <- setdiff(MEMBER_ITEMS, names(dt))
if (length(missing_member) > 0) {
  stop("Missing membership items:\n- ", paste(missing_member, collapse = "\n- "))
}

dt <- dt %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(MEMBER_ITEMS), ~ to01(ifelse(to_na(.x) < 2, 0, 1)))) %>%
  dplyr::rename(
    nhg          = c12_01,
    religious    = c12_02,
    sport        = c12_03,
    charity      = c12_04,
    political    = c12_05,
    professional = c12_06,
    union        = c12_07,
    student      = c12_08
  )

membership_vars <- c("nhg","religious","political","union","professional","charity","sport","student")

# Ensure 0/1
dt <- dt %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(membership_vars), ~ to01(.x)))

# ------------------------------------------------------------------------------
# 4) Trust outcomes (kept for dt_states; not used in LMM itself)
# ------------------------------------------------------------------------------
if (!"c02" %in% names(dt)) stop("Missing c02 (generalized trust item).")

dt <- dt %>%
  dplyr::mutate(
    trust = dplyr::case_when(
      to_na(c02) == 1 ~ 1L,
      to_na(c02) == 2 ~ 0L,
      to_na(c02) == 3 ~ 1L,
      TRUE ~ NA_integer_
    )
  )

# Neighborhood trust: choose explicitly in 00_setup.R
if (!exists("TRUST_NH_VAR")) TRUST_NH_VAR <- "c03"
if (TRUST_NH_VAR %in% names(dt)) {
  dt <- dt %>% dplyr::mutate(trust_nh = to01(to_na(.data[[TRUST_NH_VAR]])))
} else {
  dt <- dt %>% dplyr::mutate(trust_nh = NA_integer_)
  warning("Neighborhood trust var not found: ", TRUST_NH_VAR, ". Set TRUST_NH_VAR in 00_setup.R.")
}

# ------------------------------------------------------------------------------
# 5) Minimal covariates (current pipeline)
# NOTE: employment is likely not present under that name in ELSOC; it will be NA.
#       You can map it later from the correct raw variable name.
# ------------------------------------------------------------------------------
dt <- dt %>%
  dplyr::mutate(
    woman = dplyr::case_when(m0_sexo == 1 ~ 0L, m0_sexo == 2 ~ 1L, TRUE ~ NA_integer_),
    edad  = to_na(m0_edad),
    education = to_na(m01),
    employment = if ("employment" %in% names(dt)) to01(to_na(employment)) else NA_integer_
  ) %>%
  dplyr::select(id, ola, dplyr::all_of(membership_vars), trust, trust_nh, edad, woman, education, employment)

# Convert to data.table for balanced-panel operations and dt_states construction
dt <- data.table::as.data.table(dt)

# ------------------------------------------------------------------------------
# 6) Balanced panel (exactly 3 waves per id)
# ------------------------------------------------------------------------------
data.table::setorder(dt, id, ola)
panel_n <- dt[, .N, by = id]
dt <- dt[id %in% panel_n[N == 3, id]]

# ------------------------------------------------------------------------------
# 7) Build Y array (subjects x items x occasions) then convert to S (n x TT x r)
# ------------------------------------------------------------------------------
ids   <- unique(dt$id)
waves <- sort(unique(dt$ola))
subjects  <- length(ids)
items     <- length(membership_vars)
occasions <- length(waves)

if (occasions != 3) warning("Expected 3 occasions; got ", occasions, ". Check WAVES_RAW mapping.")

Y <- array(NA_integer_, dim = c(subjects, items, occasions))

for (t in seq_along(waves)) {
  dt_t <- dt[ola == waves[t]]
  dt_t <- dt_t[match(ids, id)]
  Y[, , t] <- as.matrix(dt_t[, ..membership_vars])
}

# LMest expects S as n x TT x r
S  <- aperm(Y, c(1, 3, 2))
yv <- rep(1, dim(S)[1])

# Guardrail: categories must be 0/1 or NA
stopifnot(all(is.na(S) | S %in% c(0L, 1L)))

# ------------------------------------------------------------------------------
# 8) Fit K=1..5 with MULTI-START and save fit table
# ------------------------------------------------------------------------------
fit_list <- vector("list", max(K_GRID))
fit_stats <- tibble::tibble(K = integer(), logLik = double(), npar = integer(), AIC = double(), BIC = double())

for (K in K_GRID) {
  set.seed(123)  # reproducible sequence of random starts
  fit <- fit_best_lmest(S, yv, k = K, nstart = N_STARTS, mod = MOD_TRANS, tol = TOL, maxit = MAXIT)
  fit_list[[K]] <- fit

  fit_stats <- dplyr::bind_rows(
    fit_stats,
    tibble::tibble(K = K, logLik = fit$lk, npar = fit$np, AIC = fit$aic, BIC = fit$bic)
  )

  message("Finished K=", K, " | logLik=", round(fit$lk, 2), " | BIC=", round(fit$bic, 2))
}

dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(fit_stats, here::here("output", "fit_table.csv"))

# ------------------------------------------------------------------------------
# 9) Select baseline model (K=3 unless you change K_BASELINE)
# ------------------------------------------------------------------------------
mK <- fit_list[[K_BASELINE]]
if (is.null(mK)) stop("Baseline model K=", K_BASELINE, " was not fitted.")

# ------------------------------------------------------------------------------
# 10) Conditional response probabilities: P(Y=1 | state) — ROBUST extraction
# ------------------------------------------------------------------------------
psi <- mK$Psi
if (is.null(psi)) stop("mK$Psi not found. Cannot compute profile probabilities.")

d <- dim(psi)
k <- K_BASELINE
r <- length(membership_vars)
m <- 2  # binary categories

dim_cat   <- which(d == m)[1]
dim_items <- which(d == r)[1]
dim_state <- which(d == k)[1]

if (any(is.na(c(dim_cat, dim_items, dim_state)))) {
  stop(
    "Could not identify Psi dimensions. dim(Psi) = ", paste(d, collapse=" x "),
    " | expected one dim=2 (categories), one dim=", r, " (items), one dim=", k, " (states)."
  )
}

psi_std <- aperm(psi, c(dim_items, dim_state, dim_cat))
prob1 <- psi_std[, , 2, drop = TRUE]  # items x states

prob_df <- as.data.frame(prob1)
colnames(prob_df) <- paste0("State_", 1:k)
prob_df$item <- membership_vars

prob_long <- prob_df |>
  tidyr::pivot_longer(dplyr::starts_with("State_"), names_to = "state", values_to = "prob") |>
  dplyr::mutate(
    state = dplyr::recode(state,
      "State_1" = "α (isolation)",
      "State_2" = "β (clustering)",
      "State_3" = "γ (bridging)"
    ),
    item = factor(item, levels = c("nhg","religious","political","union","professional","charity","sport","student"))
  )

p_profiles_dot <- ggplot2::ggplot(prob_long, ggplot2::aes(x = item, y = prob, color = state)) +
  ggplot2::geom_point(size = 4) +
  ggplot2::geom_segment(
    ggplot2::aes(x = item, xend = item, y = 0, yend = prob),
    linewidth = 0.6, alpha = 0.35
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = paste0("Latent portfolio profiles (K=", K_BASELINE, ", mod=", MOD_TRANS, ")"),
    x = "Association domain",
    y = "Pr(membership)"
  ) +
  theme_ssr_big() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

ggsave_safe("plot_latentclass.png", p_profiles_dot, width = 10, height = 6)
ggplot2::ggsave(here::here("output","plot_latentclass.pdf"), p_profiles_dot, width = 10, height = 6)

readr::write_csv(
  prob_df,
  here::here("output", paste0("latent_profiles_K", k, "_mod", MOD_TRANS, ".csv"))
)

p_profiles_facet <- ggplot2::ggplot(prob_long, ggplot2::aes(x = item, y = prob)) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::facet_wrap(~state, nrow = 1) +
  ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = paste0("Latent portfolio profiles (K=", K_BASELINE, ", mod=", MOD_TRANS, ")"),
    x = "Association domain",
    y = "Pr(membership)"
  ) +
  theme_ssr_big() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
                 legend.position = "none")

ggsave_safe("plot_latentclass_facet.png", p_profiles_facet, width = 12, height = 5.2)
ggplot2::ggsave(here::here("output","plot_latentclass_facet.pdf"), p_profiles_facet, width = 12, height = 5.2)

# ------------------------------------------------------------------------------
# 11) Transition matrix (average)
# ------------------------------------------------------------------------------
Pi <- mK$Pi
if (is.null(Pi)) stop("mK$Pi not found. Cannot compute transition matrix.")

Pi_avg <- apply(mK$Pi, c(1, 2), mean)

Pi_heat <- as.data.frame(Pi_avg) |>
  dplyr::mutate(from = factor(1:K_BASELINE)) |>
  tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "p") |>
  dplyr::mutate(
    to = factor(gsub("V", "", to)),
    from = dplyr::recode(from, "1" = "α", "2" = "β", "3" = "γ"),
    to   = dplyr::recode(to,   "1" = "α", "2" = "β", "3" = "γ"),
    lab = sprintf("%.2f", p)
  )

p_trans <- ggplot2::ggplot(Pi_heat, ggplot2::aes(x = to, y = from, fill = p)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = lab), size = 3.8) +
  ggplot2::scale_fill_gradient(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = paste0("Transition probabilities (average) (K = ", K_BASELINE, ", mod=", MOD_TRANS, ")"),
    x = "To state",
    y = "From state"
  ) +
  theme_ssr() +
  ggplot2::theme(legend.position = "right")

ggsave_safe("plot_transition.png", p_trans, width = 6.6, height = 5.4)
ggplot2::ggsave(here::here("output","plot_transition.pdf"), p_trans, width = 6.6, height = 5.4)

readr::write_csv(
  as.data.frame(Pi_avg) %>% dplyr::mutate(from_state = 1:k) %>% dplyr::relocate(from_state),
  here::here("output", paste0("transition_matrix_avg_K", k, "_mod", MOD_TRANS, ".csv"))
)

# ------------------------------------------------------------------------------
# 12) Most likely state per wave and α/β/γ labeling (ROBUST)
# ------------------------------------------------------------------------------
V <- mK$V
if (is.null(V)) stop("mK$V not found. Cannot assign most likely states.")

n  <- length(ids)
k  <- K_BASELINE
TT <- length(waves)

V_std <- reshape_V_to_nkT(V, n = n, k = k, TT = TT)

# Modal state per wave (Uhat: n x TT) — NOTE drop=TRUE => matrix
Uhat <- sapply(seq_len(TT), function(t) {
  max.col(V_std[, , t, drop = TRUE], ties.method = "first")
})
Uhat <- as.matrix(Uhat)  # n x TT

Uhat_dt <- data.table::data.table(
  id = ids,
  state_t1 = Uhat[, 1],
  state_t2 = Uhat[, 2],
  state_t3 = Uhat[, 3]
)

# Build dt_states
dt_states <- dt[, .(id, ola, trust, trust_nh, edad, woman, education, employment,
                    nhg, religious, political, union, professional, charity, sport, student)]
dt_states[, state := NA_integer_]

dt_states[ola == waves[1], state := Uhat_dt$state_t1[match(id, Uhat_dt$id)]]
dt_states[ola == waves[2], state := Uhat_dt$state_t2[match(id, Uhat_dt$id)]]
dt_states[ola == waves[3], state := Uhat_dt$state_t3[match(id, Uhat_dt$id)]]

dt_states[, position := factor(state, levels = 1:k, labels = c("alpha", "beta", "gamma"))]

# Sanity checks: should NOT be all alpha
pos_counts <- dt_states[, .N, by = position][order(-N)]
print(pos_counts)

if (nrow(pos_counts) < 2) {
  stop("Degenerate classification: position has <2 levels. Check V reshaping / model fit.")
}

saveRDS(dt_states, here::here("data", "dt_states.rds"))


# ------------------------------------------------------------------------------
# 13) Posterior probabilities (classification uncertainty) — SAVE STANDARDIZED
# ------------------------------------------------------------------------------
# We save:
#  (a) V_std: n x k x TT (respondent x state x wave-index)
#  (b) post_long: id-ola with p_alpha/p_beta/p_gamma

post_dir <- here::here("data")
dir.create(post_dir, showWarnings = FALSE, recursive = TRUE)

# Standardized posterior array
post_std <- V_std  # n x k x TT (already standardized above)

# Build long posteriors (id x ola)
post_long <- rbindlist(lapply(seq_len(TT), function(t){
  out <- data.table(
    id  = ids,
    ola = waves[t],
    p_alpha = post_std[, 1, t],
    p_beta  = post_std[, 2, t],
    p_gamma = post_std[, 3, t]
  )
  out
}))

# Merge posteriors into dt_states (so downstream scripts are single-input)
dt_states <- merge(dt_states, post_long, by = c("id","ola"), all.x = TRUE)

# Guardrail
if (anyNA(dt_states$p_gamma)) warning("Some p_gamma are NA after merge. Check id/ola alignment.")

saveRDS(dt_states, here::here("data", "dt_states.rds"))
saveRDS(list(ids = ids, waves = waves, V_std = post_std, post_long = post_long),
        here::here("data", "posterior_probs_std.rds"))

message("Saved dt_states with membership vars + posteriors: data/dt_states.rds")
message("Saved standardized posteriors: data/posterior_probs_std.rds")