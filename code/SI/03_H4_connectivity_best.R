# ==============================================================================
# 03_H4_connectivity_figures_themeSSR.R — H4 figures (ROBUST + BOOTSTRAP CI)
# ------------------------------------------------------------------------------
# No weights. X axis shows only existing years (2016, 2018, 2022).
# Adds bootstrap 95% CIs (cluster bootstrap by individual id) for:
#   - prev_gamma_soft   (mean posterior p_gamma)
#   - prev_gamma_strict (share gamma among strict assignments)
#   - p_two_plus        (share with >=2 domains)
#   - mean_jacc         (mean pairwise Jaccard across domain pairs)
#
# Outputs:
#   - output/H4_connectivity_by_wave_jaccard.csv
#   - output/H4_connectivity_by_wave_bootstrap.csv
#   - output/H4_panel.png/.pdf + component plots
#   - output/H4_connectivity_log.txt
#
# Requires:
#   code/00_setup.R
#   data/dt_states_cov_FIXED.rds OR data/dt_states_cov.rds OR data/dt_states.rds
#   data/posterior_probs_cov_std.rds (only if posteriors need rebuild)
#   data/dt_analysis.rds (only if domain vars missing in dt_states)
# ==============================================================================

source(here::here("code","00_setup.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(readr)
  library(patchwork)
})

dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)
dir.create(here::here("data"),   showWarnings = FALSE, recursive = TRUE)

LOG_PATH <- here::here("output", "H4_connectivity_log.txt")
file.remove(LOG_PATH)  # reset log each run

log_line <- function(...) {
  txt <- paste0(..., collapse = "")
  cat(txt, "\n")
  write(txt, file = LOG_PATH, append = TRUE)
}

# ------------------------------------------------------------------------------
# 1) Load dt_states (prefer FIXED > cov > nocov)
# ------------------------------------------------------------------------------
path_fixed <- here::here("data","dt_states_cov_FIXED.rds")
path_cov   <- here::here("data","dt_states_cov.rds")
path_noc   <- here::here("data","dt_states.rds")

if (file.exists(path_fixed)) {
  dt <- readRDS(path_fixed)
  log_line("Loaded: ", path_fixed)
} else if (file.exists(path_cov)) {
  dt <- readRDS(path_cov)
  log_line("Loaded: ", path_cov)
} else if (file.exists(path_noc)) {
  dt <- readRDS(path_noc)
  log_line("Loaded: ", path_noc)
} else {
  stop("No dt_states file found in data/. Expected one of: dt_states_cov_FIXED.rds, dt_states_cov.rds, dt_states.rds")
}
dt <- as.data.table(dt)

# ------------------------------------------------------------------------------
# 2) Standardize keys (id as character; ola numeric)
# ------------------------------------------------------------------------------
if (!("id" %in% names(dt)) && ("idencuesta" %in% names(dt))) setnames(dt, "idencuesta", "id")
if (!("ola" %in% names(dt))) stop("Missing 'ola' column in dt_states.")
dt[, id := as.character(id)]
dt[, ola := as.integer(ola)]
dt[, wave := ola]

if ("position" %in% names(dt))        dt[, position := as.character(position)]
if ("position_strict" %in% names(dt)) dt[, position_strict := as.character(position_strict)]

# ------------------------------------------------------------------------------
# 3) Ensure domain vars exist (merge from dt_analysis if needed)
# ------------------------------------------------------------------------------
domains_named <- unique(c(DOMAIN_VARS, "otra"))
domains_c12   <- paste0("c12_", sprintf("%02d", 1:9))

has_named <- all(DOMAIN_VARS %in% names(dt))
has_c12   <- any(domains_c12 %in% names(dt))

if (!has_named && !has_c12) {
  anal_path <- here::here("data","dt_analysis.rds")
  if (!file.exists(anal_path)) stop("Domain vars not found in dt_states and dt_analysis.rds not available.")
  da <- as.data.table(readRDS(anal_path))

  if (!("id" %in% names(da)) && ("idencuesta" %in% names(da))) setnames(da, "idencuesta", "id")
  if (!("id" %in% names(da)) || !("ola" %in% names(da))) stop("dt_analysis.rds missing id/ola keys.")

  da[, id := as.character(id)]
  da[, ola := as.integer(ola)]

  keep_dom <- intersect(names(da), domains_named)
  if (length(keep_dom) == 0) stop("dt_analysis.rds does not contain named domain vars.")

  setkey(dt, id, ola)
  setkey(da, id, ola)

  # Add domains to dt (preserve dt rows)
  dt <- da[dt]
  log_line("Merged domain vars from dt_analysis.rds: ", paste(keep_dom, collapse=", "))
}

# Decide domain set
if (all(DOMAIN_VARS %in% names(dt))) {
  domains <- DOMAIN_VARS
  if ("otra" %in% names(dt)) domains <- c(domains, "otra")
  log_line("Using named domains: ", paste(domains, collapse=", "))
} else if (any(domains_c12 %in% names(dt))) {
  domains <- domains_c12[domains_c12 %in% names(dt)]
  log_line("Using c12_* domains: ", paste(domains, collapse=", "))
} else {
  stop("Could not detect domain variables after merge attempts.")
}

# ------------------------------------------------------------------------------
# 4) Posterior validation
# ------------------------------------------------------------------------------
posteriors_ok <- function(DT, tol = 1e-6) {
  if (!all(c("p_alpha","p_beta","p_gamma") %in% names(DT))) return(FALSE)

  DT[, p_alpha := as.numeric(p_alpha)]
  DT[, p_beta  := as.numeric(p_beta)]
  DT[, p_gamma := as.numeric(p_gamma)]

  bad_range <- DT[, any(p_alpha < -tol | p_beta < -tol | p_gamma < -tol |
                        p_alpha > 1+tol | p_beta > 1+tol | p_gamma > 1+tol,
                        na.rm = TRUE)]
  if (isTRUE(bad_range)) return(FALSE)

  md <- DT[, max(abs((p_alpha + p_beta + p_gamma) - 1), na.rm = TRUE)]
  is.finite(md) && md <= tol
}

# ------------------------------------------------------------------------------
# 5) Rebuild posteriors if invalid/missing (and SAVE FIXED)
# ------------------------------------------------------------------------------
if (!posteriors_ok(dt)) {
  log_line("Posteriors missing/invalid in dt_states. Rebuilding from V_std...")

  obj_path <- here::here("data","posterior_probs_cov_std.rds")
  if (!file.exists(obj_path)) stop("Cannot rebuild posteriors: missing data/posterior_probs_cov_std.rds")
  obj <- readRDS(obj_path)

  V_std <- obj$V_std
  ids   <- as.character(obj$ids)
  waves <- as.integer(obj$waves)

  if (!is.null(obj$mapping)) {
    alpha_state <- obj$mapping$alpha
    beta_state  <- obj$mapping$beta
    gamma_state <- obj$mapping$gamma
  } else {
    alpha_state <- obj$alpha_state
    beta_state  <- obj$beta_state
    gamma_state <- obj$gamma_state
  }
  if (any(sapply(list(alpha_state,beta_state,gamma_state), is.null))) {
    stop("Missing alpha/beta/gamma mapping in posterior_probs_cov_std.rds")
  }

  THRESH_STRICT <- obj$THRESH_STRICT
  if (is.null(THRESH_STRICT) || !is.finite(THRESH_STRICT)) THRESH_STRICT <- 0.60

  drop_cols <- intersect(
    names(dt),
    c("p_alpha","p_beta","p_gamma","p_max","state_num","position","position_strict")
  )
  if (length(drop_cols) > 0) dt[, (drop_cols) := NULL]

  dt[, key := paste(id, ola)]
  if (dt[, any(duplicated(key))]) {
    log_line("WARNING: duplicate id-ola rows; keeping first.")
    dt <- dt[!duplicated(key)]
  }
  dt[, key := NULL]

  post_long <- rbindlist(lapply(seq_along(waves), function(t) {
    pp <- V_std[, , t, drop = TRUE]
    data.table(
      id = ids,
      ola = waves[t],
      p_alpha = as.numeric(pp[, alpha_state]),
      p_beta  = as.numeric(pp[, beta_state]),
      p_gamma = as.numeric(pp[, gamma_state]),
      p_max   = as.numeric(apply(pp, 1, max))
    )
  }))

  modal_long <- rbindlist(lapply(seq_along(waves), function(t) {
    pp <- V_std[, , t, drop = TRUE]
    st <- max.col(pp, ties.method = "first")
    pos <- ifelse(st == alpha_state, "alpha",
                  ifelse(st == gamma_state, "gamma", "beta"))
    data.table(
      id = ids,
      ola = waves[t],
      state_num = as.integer(st),
      position  = as.character(pos)
    )
  }))

  setkey(dt, id, ola)
  setkey(post_long, id, ola)
  setkey(modal_long, id, ola)

  dt <- post_long[dt]    # preserve dt rows
  dt <- modal_long[dt]

  dt[, position_strict := position]
  dt[p_max < THRESH_STRICT, position_strict := NA_character_]
  dt[, wave := ola]

  md <- dt[, max(abs((p_alpha + p_beta + p_gamma) - 1), na.rm = TRUE)]
  log_line("Rebuild check max_dev=", signif(md, 4), " | THRESH_STRICT=", THRESH_STRICT)
  if (!is.finite(md) || md > 1e-6) stop("Rebuild failed: posteriors still invalid.")

  saveRDS(as.data.frame(dt), here::here("data","dt_states_cov_FIXED.rds"))
  log_line("Saved FIXED: data/dt_states_cov_FIXED.rds")
}

if (!("p_max" %in% names(dt))) dt[, p_max := pmax(p_alpha, p_beta, p_gamma, na.rm = TRUE)]
if (!("position_strict" %in% names(dt))) {
  THRESH_STRICT <- 0.60
  dt[, position_strict := as.character(position)]
  dt[p_max < THRESH_STRICT, position_strict := NA_character_]
}

# ------------------------------------------------------------------------------
# 6) Map wave -> year (and keep only existing years on x axis)
# ------------------------------------------------------------------------------
# Your pipeline uses ola recoded to 1,2,3 OR raw 1,3,6. Support both.
dt[, year := fifelse(ola %in% c(1L,2L,3L), WAVES_YEARS[ola],
                     fifelse(ola == 1L, WAVES_YEARS[1],
                             fifelse(ola == 3L, WAVES_YEARS[2],
                                     fifelse(ola == 6L, WAVES_YEARS[3], NA_integer_))))]

years_present <- sort(unique(dt$year[!is.na(dt$year)]))
if (length(years_present) < 2) stop("Could not map waves to years. Check ola coding.")
log_line("Years present: ", paste(years_present, collapse=", "))

# ------------------------------------------------------------------------------
# 7) H4 metric computation
# ------------------------------------------------------------------------------
mean_pairwise_jaccard <- function(X) {
  P <- crossprod(X)
  diagv <- diag(P)
  d <- ncol(X)

  Jvals <- numeric(0)
  for (a in 1:(d-1)) {
    for (b in (a+1):d) {
      union_ab <- diagv[a] + diagv[b] - P[a,b]
      if (union_ab > 0) Jvals <- c(Jvals, P[a,b] / union_ab)
    }
  }
  if (length(Jvals) == 0) return(NA_real_)
  mean(Jvals, na.rm = TRUE)
}

compute_h4_by_year <- function(DT) {
  DT[, {
    X <- as.matrix(.SD)
    X[is.na(X)] <- 0
    storage.mode(X) <- "numeric"

    dom_count  <- rowSums(X)
    .(
      mean_jacc = mean_pairwise_jaccard(X),
      p_two_plus = mean(dom_count >= 2),

      prev_gamma_soft   = mean(p_gamma, na.rm = TRUE),
      prev_gamma_hard   = mean(position == "gamma", na.rm = TRUE),
      prev_gamma_strict = mean(position_strict == "gamma", na.rm = TRUE),

      mean_pmax   = mean(p_max, na.rm = TRUE),
      share_strict = mean(!is.na(position_strict))
    )
  }, by = year, .SDcols = domains][order(year)]
}

h4_point <- compute_h4_by_year(dt)
write_csv(as.data.frame(h4_point), here::here("output","H4_connectivity_by_wave_jaccard.csv"))
log_line("Saved: output/H4_connectivity_by_wave_jaccard.csv")

# ------------------------------------------------------------------------------
# 8) Bootstrap 95% CI (cluster bootstrap by id, within year)
# ------------------------------------------------------------------------------
set.seed(123)
B <- 500  # increase to 1000 if you want smoother intervals

boot_one_year <- function(DT_year) {
  # cluster bootstrap: resample ids with replacement
  ids <- unique(DT_year$id)
  draw <- sample(ids, size = length(ids), replace = TRUE)

  # build bootstrap sample by stacking clusters; keep multiplicity
  bs <- rbindlist(lapply(draw, function(i) DT_year[id == i]), use.names = TRUE, fill = TRUE)
  # compute metrics on bootstrap sample
  compute_h4_by_year(bs)[, year := unique(DT_year$year)][1]
}

boot_all <- rbindlist(lapply(years_present, function(yy) {
  DTy <- dt[year == yy]
  reps <- rbindlist(lapply(seq_len(B), function(b) boot_one_year(DTy)), fill = TRUE)
  reps[, .(
    year = yy,
    prev_gamma_soft_lo   = quantile(prev_gamma_soft,   0.025, na.rm=TRUE),
    prev_gamma_soft_hi   = quantile(prev_gamma_soft,   0.975, na.rm=TRUE),
    prev_gamma_strict_lo = quantile(prev_gamma_strict, 0.025, na.rm=TRUE),
    prev_gamma_strict_hi = quantile(prev_gamma_strict, 0.975, na.rm=TRUE),
    p_two_plus_lo        = quantile(p_two_plus,        0.025, na.rm=TRUE),
    p_two_plus_hi        = quantile(p_two_plus,        0.975, na.rm=TRUE),
    mean_jacc_lo         = quantile(mean_jacc,         0.025, na.rm=TRUE),
    mean_jacc_hi         = quantile(mean_jacc,         0.975, na.rm=TRUE)
  )]
}), fill = TRUE)

h4 <- merge(h4_point, boot_all, by = "year", all.x = TRUE)
setorder(h4, year)

write_csv(as.data.frame(h4), here::here("output","H4_connectivity_by_wave_bootstrap.csv"))
log_line("Saved: output/H4_connectivity_by_wave_bootstrap.csv (B=", B, ")")

# ------------------------------------------------------------------------------
# 9) Plot helpers: fixed x ticks and honest y scaling
# ------------------------------------------------------------------------------
THEME <- theme_ssr()

x_scale_years <- scale_x_continuous(breaks = years_present, labels = years_present)

scale_pct0 <- function() {
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.06))
  )
}

# For Jaccard: keep within [0, NA] but do not force to 0% scale; still start at 0 for interpretability
scale_jacc0 <- function() {
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.06))
  )
}

# ------------------------------------------------------------------------------
# 10) Figures (with 95% CI error bars)
# ------------------------------------------------------------------------------
p_gamma_soft <- ggplot(h4, aes(x = year, y = prev_gamma_soft, group = 1)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = prev_gamma_soft_lo, ymax = prev_gamma_soft_hi), width = 0.15) +
  x_scale_years +
  scale_pct0() +
  labs(x="Wave", y="Mean posterior Pr(γ)", title="Bridging prevalence by wave (soft)") +
  THEME

p_twoplus <- ggplot(h4, aes(x = year, y = p_two_plus, group = 1)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = p_two_plus_lo, ymax = p_two_plus_hi), width = 0.15) +
  x_scale_years +
  scale_pct0() +
  labs(x="Wave", y="Share with ≥2 domains", title="Multi-domain participation (≥2 domains)") +
  THEME

p_gamma_strict <- ggplot(h4, aes(x = year, y = prev_gamma_strict, group = 1)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = prev_gamma_strict_lo, ymax = prev_gamma_strict_hi), width = 0.15) +
  x_scale_years +
  scale_pct0() +
  labs(x="Wave", y="Pr(γ) among strict assignments", title="Bridging prevalence by wave (strict)") +
  THEME

p_jacc <- ggplot(h4, aes(x = year, y = mean_jacc, group = 1)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = mean_jacc_lo, ymax = mean_jacc_hi), width = 0.15) +
  x_scale_years +
  scale_jacc0() +
  labs(x="Wave", y="Mean pairwise Jaccard (domain overlap)", title="Portfolio-implied cross-domain connectivity") +
  THEME

p_panel <- (p_gamma_soft | p_twoplus) / (p_gamma_strict | p_jacc) +
  plot_annotation(
    title = "H4: Portfolio dynamics and portfolio-implied connectivity (with bootstrap 95% CI)",
    theme = THEME
  )

# ------------------------------------------------------------------------------
# 11) Save
# ------------------------------------------------------------------------------
ggsave_safe("H4_prev_gamma_soft.png",   p_gamma_soft,   width = 7.2, height = 4.6)
ggsave_safe("H4_p_two_plus.png",        p_twoplus,      width = 7.2, height = 4.6)
ggsave_safe("H4_prev_gamma_strict.png", p_gamma_strict, width = 7.2, height = 4.6)
ggsave_safe("H4_mean_jacc.png",         p_jacc,         width = 7.2, height = 4.6)
ggsave_safe("H4_panel.png",             p_panel,        width = 10.8, height = 8.2)

ggsave(here::here("output","H4_prev_gamma_soft.pdf"),   p_gamma_soft,   width = 7.2, height = 4.6)
ggsave(here::here("output","H4_p_two_plus.pdf"),        p_twoplus,      width = 7.2, height = 4.6)
ggsave(here::here("output","H4_prev_gamma_strict.pdf"), p_gamma_strict, width = 7.2, height = 4.6)
ggsave(here::here("output","H4_mean_jacc.pdf"),         p_jacc,         width = 7.2, height = 4.6)
ggsave(here::here("output","H4_panel.pdf"),             p_panel,        width = 10.8, height = 8.2)

log_line("DONE: H4 figures + bootstrap CI saved in output/")
log_line("Log: ", LOG_PATH)