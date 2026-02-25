# ==============================================================================
# 02_latent_markov.R — Latent Markov model for membership portfolios (α/β/γ)
# ------------------------------------------------------------------------------
# Inputs:
#   data/ELSOC_Long.RData  (expects object ELSOC_Long)
# Outputs (output/):
#   - fit_table.csv
#   - plot_latentclass.png
#   - plot_transition.png
#   - dt_states.rds (long panel with positions per wave)
#   - posterior_probs.rds (if available from LMest)
# ==============================================================================

source(here::here("code", "00_setup.R"))

# ---- Load data ----
stop_if_missing(c(here::here("data", "ELSOC_Long.RData")))
load(here::here("data", "ELSOC_Long.RData"))  # loads ELSOC_Long

dt <- data.table::as.data.table(ELSOC_Long)

# ---- Keep only relevant columns (adapt here if your variable names differ) ----
keep <- c(
  "id", "ola",
  "nhg","religious","political","union","professional","charity","sport","student",
  "trust","trust_nh",
  "edad","woman","education","employment"
)
missing_cols <- setdiff(keep, names(dt))
if (length(missing_cols) > 0) {
  stop("ELSOC_Long is missing expected columns:\n- ", paste(missing_cols, collapse = "\n- "))
}
dt <- dt[, ..keep]

# ---- Recode membership indicators to 0/1 ----
membership_vars <- c("nhg","religious","political","union","professional","charity","sport","student")
for (v in membership_vars) dt[, (v) := to01(get(v))]

# Trust outcomes to 0/1 (conservative: treat only explicit trust as 1)
dt[, trust    := to01(trust)]
dt[, trust_nh := to01(trust_nh)]

# ---- Balanced panel (3 waves) ----
data.table::setorder(dt, id, ola)
panel_n <- dt[, .N, by = id]
dt <- dt[id %in% panel_n[N == 3, id]]

# ---- Build Y array for LMest: (subjects, items, occasions) ----
ids <- unique(dt$id)
waves <- sort(unique(dt$ola))
subjects  <- length(ids)
items     <- length(membership_vars)
occasions <- length(waves)

Y <- array(NA_integer_, dim = c(subjects, items, occasions))

for (t in seq_along(waves)) {
  dt_t <- dt[ola == waves[t]]
  dt_t <- dt_t[match(ids, id)]
  Y[, , t] <- as.matrix(dt_t[, ..membership_vars])
}

# ---- Fit K=1..5 and select by BIC ----
fit_list <- vector("list", 5)
fit_stats <- tibble::tibble(K = integer(), logLik = double(), npar = integer(), AIC = double(), BIC = double())

for (K in 1:5) {
  set.seed(123)
  fit <- LMest::est_lm_basic(Y = Y, k = K, out_se = FALSE)
  fit_list[[K]] <- fit
  fit_stats <- dplyr::bind_rows(
    fit_stats,
    tibble::tibble(K = K, logLik = fit$lk, npar = fit$np, AIC = fit$aic, BIC = fit$bic)
  )
}

readr::write_csv(fit_stats, here::here("output", "fit_table.csv"))

# ---- Select K=3 (paper baseline) ----
m3 <- fit_list[[3]]

# ---- Conditional response probs: P(Y=1|state) ----
# LMest stores Psi as items x states x categories; for binary, categories=2 and Y=1 corresponds to category 2
prob1 <- m3$Psi[, , 2, drop = TRUE]
prob_df <- as.data.frame(prob1)
colnames(prob_df) <- paste0("State_", 1:3)
prob_df$item <- membership_vars
prob_df <- dplyr::relocate(prob_df, item)

prob_long <- prob_df |>
  tidyr::pivot_longer(cols = starts_with("State_"), names_to = "state", values_to = "prob")

p_profiles <- ggplot2::ggplot(prob_long, ggplot2::aes(x = item, y = prob, group = state)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    x = "Association domain",
    y = "P(membership=1 | latent state)",
    title = "Latent portfolio profiles (K = 3)"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggsave_safe("plot_latentclass.png", p_profiles, width = 9, height = 5)

# ---- Transition matrix (average across transitions) ----
Pi <- m3$Pi  # states x states x (T-1)
Pi_avg <- apply(Pi, c(1,2), mean)

Pi_long <- as.data.frame(Pi_avg) |>
  dplyr::mutate(from_state = factor(1:3)) |>
  tidyr::pivot_longer(cols = -from_state, names_to = "to_state", values_to = "prob")

p_trans <- ggplot2::ggplot(Pi_long, ggplot2::aes(x = from_state, y = prob, fill = to_state)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Average transition probabilities",
    x = "From state",
    y = "Probability",
    fill = "To state"
  )

ggsave_safe("plot_transition.png", p_trans, width = 7, height = 5)

# ---- Assign positions α/β/γ (most likely state, wave-specific) ----
Uhat <- m3$Ul  # subjects x occasions (most likely latent state per wave)

Uhat_dt <- data.table::data.table(
  id = ids,
  state_t1 = Uhat[,1],
  state_t2 = Uhat[,2],
  state_t3 = Uhat[,3]
)

dt_states <- dt[, .(id, ola, trust, trust_nh, edad, woman, education, employment)]
dt_states[, state := NA_integer_]

dt_states[ola == waves[1], state := Uhat_dt$state_t1[match(id, Uhat_dt$id)]]
dt_states[ola == waves[2], state := Uhat_dt$state_t2[match(id, Uhat_dt$id)]]
dt_states[ola == waves[3], state := Uhat_dt$state_t3[match(id, Uhat_dt$id)]]

# IMPORTANT: mapping of latent states to α/β/γ should follow your empirical interpretation.
# By default we assume: 1=α (isolation), 2=β (clustering), 3=γ (bridging).
dt_states[, position := factor(
  state,
  levels = 1:3,
  labels = c("alpha", "beta", "gamma")
)]

# Save main analysis panel with positions
saveRDS(dt_states, here::here("data", "dt_states.rds"))

# ---- Save posterior probabilities if available (classification uncertainty) ----
# Depending on LMest version, posteriors may be in m3$V (subjects x states x occasions).
post <- NULL
if (!is.null(m3$V)) {
  post <- m3$V
}
if (!is.null(post)) {
  saveRDS(list(ids = ids, waves = waves, post = post), here::here("data", "posterior_probs.rds"))
}
