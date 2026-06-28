# =============================================================================
# 03g_mediation_analysis.R
# Mediation: Demographics → Portfolio Class → Trust
# Method: product-of-coefficients + delta method + parametric bootstrap (B=5000)
#
# a-path: AME of demographic X on P(class = k), estimated via OLS with
#         cluster-robust SEs (clustered by individual) on LMM posteriors
# b-path: AME of class k vs γ on trust, from stored Stage III glmer objects
# IE_k  = a_k × b_k  (for k = β, α)
# Total IE = IE_β + IE_α
#
# Focal demographics: edu_bin (0=low, 1=high), woman (0=man, 1=woman)
# Age enters as control in a-path regressions, not as focal mediator.
# Reference class: γ (Bridging), clase = 2
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(marginaleffects)
  library(sandwich)
  library(lmtest)
})

set.seed(20240628)
B <- 5000L   # parametric bootstrap draws

# ── 0. Paths ──────────────────────────────────────────────────────────────────
DATA_DIR   <- "data"
OUTPUT_DIR <- "output"

# ── 1. Load data and models ───────────────────────────────────────────────────
mobj  <- readRDS(file.path(DATA_DIR, "trust_model_objects.rds"))
dt_sc <- readRDS(file.path(DATA_DIR, "dt_states_cov.rds"))

# mobj$dt already contains p_alpha, p_beta, p_gamma, edad — no join needed
dt <- mobj$dt |>
  mutate(
    age_num = case_when(
      edad == "18_24" ~ 21.0,
      edad == "25_34" ~ 29.5,
      edad == "35_44" ~ 39.5,
      edad == "45_54" ~ 49.5,
      edad == "55_64" ~ 59.5,
      edad == "65"    ~ 70.0,
      TRUE            ~ NA_real_
    ),
    age_c = as.numeric(scale(age_num))
  ) |>
  filter(!is.na(p_alpha), !is.na(edu_bin), !is.na(woman), !is.na(age_c))

cat(sprintf("Working dataset: %d person-waves, %d individuals\n",
            nrow(dt), length(unique(dt$idencuesta))))

# ── 2. B-path: AMEs from stored Stage III glmer objects ───────────────────────
# t4_m2 = trust (generalized), t4_m4 = trust_nh (neighbourhood)
# Reference: γ = clase 2; contrasts: "1 - 2" (β vs γ), "3 - 2" (α vs γ)

extract_b <- function(model, outcome_lbl) {
  ame <- avg_slopes(model, variables = "clase_f", vcov = TRUE)
  df  <- as.data.frame(ame)
  data.frame(
    outcome  = outcome_lbl,
    via      = c("beta",  "alpha"),
    contrast = c("1 - 2", "3 - 2"),
    b_ame    = as.numeric(df$estimate[match(c("1 - 2","3 - 2"), df$contrast)]),
    b_se     = as.numeric(df$std.error[match(c("1 - 2","3 - 2"), df$contrast)]),
    stringsAsFactors = FALSE
  )
}

b_path <- bind_rows(
  extract_b(mobj$t4_m2, "trust"),
  extract_b(mobj$t4_m4, "trust_nh")
)

cat("\nB-path (Stage III AMEs):\n")
print(b_path)

# ── 3. A-path: OLS on LMM posteriors, cluster-robust SEs ─────────────────────
# For focal X ∈ {edu_bin, woman}:
#   lm(p_k ~ X + other_demo + age_c)  clustered by idencuesta
# AME of binary X = regression coefficient (LPM)

get_a_ame <- function(p_outcome, focal_x, data) {
  other <- setdiff(c("edu_bin", "woman", "age_c"), focal_x)
  fml   <- as.formula(
    paste(p_outcome, "~", paste(c(focal_x, other), collapse = " + "))
  )
  m   <- lm(fml, data = data)
  ct  <- coeftest(m, vcov = vcovCL(m, cluster = ~ idencuesta))
  idx <- which(rownames(ct) == focal_x)
  data.frame(
    focal_x  = focal_x,
    via      = sub("p_", "", p_outcome),   # "beta" or "alpha"
    a_ame    = ct[idx, "Estimate"],
    a_se     = ct[idx, "Std. Error"],
    stringsAsFactors = FALSE
  )
}

a_path <- bind_rows(lapply(c("edu_bin", "woman"), function(fx) {
  bind_rows(lapply(c("p_beta", "p_alpha"), function(pk) get_a_ame(pk, fx, dt)))
}))

cat("\nA-path (demographics → class posteriors):\n")
print(a_path)

# ── 4. Indirect effects: product-of-coefficients + delta method ───────────────

compute_indirect <- function(fx, out_lbl) {
  rows <- list()

  for (k in c("beta", "alpha")) {
    a_row <- a_path |> filter(focal_x == fx, via == k)
    b_row <- b_path |> filter(outcome == out_lbl, via == k)

    a <- a_row$a_ame;  se_a <- a_row$a_se
    b <- b_row$b_ame;  se_b <- b_row$b_se

    ie    <- a * b
    se_ie <- sqrt(b^2 * se_a^2 + a^2 * se_b^2)
    z     <- ie / se_ie
    pval  <- 2 * pnorm(-abs(z))

    rows[[k]] <- list(
      focal_x = fx, outcome = out_lbl, via = k,
      a_ame = a, a_se = se_a,
      b_ame = b, b_se = se_b,
      ie = ie, ie_se = se_ie, z = z, p_val = pval
    )
  }

  # Total indirect (sum over β and α paths, assuming independence)
  ie_tot  <- rows$beta$ie  + rows$alpha$ie
  se_tot  <- sqrt(rows$beta$ie_se^2 + rows$alpha$ie_se^2)
  z_tot   <- ie_tot / se_tot
  p_tot   <- 2 * pnorm(-abs(z_tot))

  rows$TOTAL <- list(
    focal_x = fx, outcome = out_lbl, via = "TOTAL",
    a_ame = NA, a_se = NA, b_ame = NA, b_se = NA,
    ie = ie_tot, ie_se = se_tot, z = z_tot, p_val = p_tot
  )

  bind_rows(lapply(rows, as.data.frame))
}

indirect <- bind_rows(
  compute_indirect("edu_bin", "trust"),
  compute_indirect("edu_bin", "trust_nh"),
  compute_indirect("woman",   "trust"),
  compute_indirect("woman",   "trust_nh")
)

# ── 5. Parametric bootstrap CIs (B = 5000) ────────────────────────────────────
# Draw a ~ N(a_ame, a_se²) and b ~ N(b_ame, b_se²) independently,
# compute IE = a × b for each draw; report 2.5th and 97.5th percentiles.

boot_ci <- function(fx, out_lbl) {
  rows <- list()

  for (k in c("beta", "alpha")) {
    a_row <- a_path |> filter(focal_x == fx, via == k)
    b_row <- b_path |> filter(outcome == out_lbl, via == k)

    a_draws <- rnorm(B, a_row$a_ame, a_row$a_se)
    b_draws <- rnorm(B, b_row$b_ame, b_row$b_se)
    ie_draws <- a_draws * b_draws

    rows[[k]] <- data.frame(
      focal_x = fx, outcome = out_lbl, via = k,
      ci_lo = quantile(ie_draws, 0.025),
      ci_hi = quantile(ie_draws, 0.975),
      stringsAsFactors = FALSE
    )
  }

  # Total
  a_b_draws <- rnorm(B, a_path$a_ame[a_path$focal_x==fx & a_path$via=="beta"],
                        a_path$a_se[ a_path$focal_x==fx & a_path$via=="beta"])
  a_a_draws <- rnorm(B, a_path$a_ame[a_path$focal_x==fx & a_path$via=="alpha"],
                        a_path$a_se[ a_path$focal_x==fx & a_path$via=="alpha"])
  b_b_draws <- rnorm(B, b_path$b_ame[b_path$outcome==out_lbl & b_path$via=="beta"],
                        b_path$b_se[ b_path$outcome==out_lbl & b_path$via=="beta"])
  b_a_draws <- rnorm(B, b_path$b_ame[b_path$outcome==out_lbl & b_path$via=="alpha"],
                        b_path$b_se[ b_path$outcome==out_lbl & b_path$via=="alpha"])

  tot_draws <- a_b_draws * b_b_draws + a_a_draws * b_a_draws
  rows$TOTAL <- data.frame(
    focal_x = fx, outcome = out_lbl, via = "TOTAL",
    ci_lo = quantile(tot_draws, 0.025),
    ci_hi = quantile(tot_draws, 0.975),
    stringsAsFactors = FALSE
  )

  bind_rows(rows)
}

boot_cis <- bind_rows(
  boot_ci("edu_bin", "trust"),
  boot_ci("edu_bin", "trust_nh"),
  boot_ci("woman",   "trust"),
  boot_ci("woman",   "trust_nh")
)

# ── 6. Combine and report ─────────────────────────────────────────────────────
sig_stars <- function(p) {
  case_when(
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    p < .10  ~ "†",
    TRUE     ~ ""
  )
}

results <- indirect |>
  left_join(boot_cis, by = c("focal_x", "outcome", "via")) |>
  mutate(stars = sig_stars(p_val))

cat("\n========================================================\n")
cat("INDIRECT EFFECTS — product-of-coefficients + bootstrap\n")
cat("========================================================\n")
print(
  results |>
    select(focal_x, outcome, via, a_ame, b_ame, ie, ie_se, ci_lo, ci_hi, p_val, stars),
  digits = 4
)

# ── 7. Save ───────────────────────────────────────────────────────────────────
write_csv(results, file.path(OUTPUT_DIR, "mediation_indirect_effects.csv"))
write_csv(a_path,  file.path(OUTPUT_DIR, "mediation_apath.csv"))
write_csv(b_path,  file.path(OUTPUT_DIR, "mediation_bpath.csv"))

cat("\nFiles saved:\n")
cat("  output/mediation_indirect_effects.csv\n")
cat("  output/mediation_apath.csv\n")
cat("  output/mediation_bpath.csv\n")
