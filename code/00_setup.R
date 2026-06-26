# ── 00_setup.R — Global setup (packages, paths, helpers) ─────────────────────

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(data.table)
})

DIR_DATA <- here::here("data")
DIR_OUT  <- here::here("output")
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# ── Global knobs ──────────────────────────────────────────────────────────────

WAVES_RAW   <- c(1, 3, 6)
WAVES_YEARS <- c(2016, 2018, 2022)

MEMBER_ITEMS <- paste0("c12_0", 1:8)
MEMBER_ITEMS_9TH <- "c12_09"   # set to NULL to exclude the 9th domain

# Membership coding: c12 == 1 (non-member), 2 (inactive), 3 (active)
# Main analysis uses any_member (c12 >= 2); active_only (c12 == 3) used in sensitivity S7
MEMBER_CODE_LOGIC <- "any_member"

# Generalized trust (c02 == 1)
TRUST_GEN_VAR  <- "c02"
TRUST_GEN_CODE <- 1L

# Neighbourhood trust: t01 >= 4 (scale 1-5)
TRUST_NH_VAR       <- "t01"
TRUST_NH_THRESHOLD <- 4L    # t01 > 3 equivalent to t01 >= 4

# Employment: m02 <= 2 (full-time or part-time employed)
EMPLOY_VAR   <- "m02"
EMPLOY_CODES <- 1L:2L

# Education recode to 5 levels (matching Stata: recode m01 1/3=1, 4=2, 5=3, 6/7=4, 8/10=5)
EDU_VAR <- "m01"
EDU_THRESHOLD_BIN <- 4L   # >= 4 in the 5-level scale = técnica or university (~34%)

# Subjective wellbeing: s01 > 4
SWB_VAR       <- "s01"
SWB_THRESHOLD <- 4L

# Couple status: m36 < 4 (living with partner)
COUPLE_VAR   <- "m36"
COUPLE_CODES <- 1L:3L

# Institutional trust (raw; control in Stata models)
TINST_VAR <- "t02_01"

# Sample filters: keep muestra == 1 and tipo_atricion == 1
FILTER_MUESTRA       <- TRUE
FILTER_TIPO_ATRICION <- TRUE

# ── Control variables ─────────────────────────────────────────────────────────
# CONTROLS_TRUST replicates estimation.do: swb, employed, couple only.
# Education is excluded from trust models to avoid collinearity with class.
CONTROLS_BASE  <- c("edad", "woman", "education")
CONTROLS_EXTRA <- c("swb", "employed", "couple", "tinst")
CONTROLS_FULL  <- c(CONTROLS_BASE, CONTROLS_EXTRA)

# ── LMM parameters ────────────────────────────────────────────────────────────
K_GRID     <- 1:5
K_BASELINE <- 3
N_STARTS   <- 30
MOD_TRANS  <- 1
TOL        <- 1e-8
MAXIT      <- 1000

N_DOMAINS <- 8L

# ── Helpers ───────────────────────────────────────────────────────────────────

stop_if_missing <- function(paths) {
  miss <- paths[!file.exists(paths)]
  if (length(miss) > 0)
    stop("Missing files:\n- ", paste(miss, collapse = "\n- "))
}

# ELSOC missing codes: -999 = no response, -888 = don't know,
#                      -777 = technical error, -666 = incomplete survey
to_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% c(-999, -888, -777, -666, 777777, 999999)] <- NA_real_
  x
}

# Two membership logics matching the ELSOC codebook
member_binary <- function(x, logic = MEMBER_CODE_LOGIC) {
  x <- to_na(x)
  if (logic == "any_member") {
    as.integer(!is.na(x) & x >= 2L)
  } else {
    as.integer(!is.na(x) & x == 3L)
  }
}
member_active <- function(x) member_binary(x)

# Employment: m02 <= 2 (full-time or part-time)
code_employment <- function(x, codes = EMPLOY_CODES) {
  x <- to_na(x)
  as.integer(!is.na(x) & x %in% codes)
}

# Generalized trust — replicates estimation.do exactly:
#   gen social_trust = 1 if c02 == 1
#   replace social_trust = 0 if c02 > 1
# In Stata, missing (.) > 1 is TRUE, so missing c02 maps to 0, not NA.
# We replicate this without calling to_na() first.
code_gen_trust <- function(x, trust_code = TRUST_GEN_CODE) {
  x_num <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    x_num == trust_code                          ~ 1L,
    !is.na(x_num) & x_num != trust_code          ~ 0L,
    is.na(x_num)                                 ~ 0L
  )
}

# Neighbourhood trust — replicates estimation.do:
#   gen trust_nb = 1 if t01 > 3; replace trust_nb = 0 if t01 <= 3
# In Stata, missing (.) is +Inf so missing t01 maps to 1.
# We match this behaviour; in practice missing t01 is rare.
code_nh_trust <- function(x, threshold = TRUST_NH_THRESHOLD) {
  x_num <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    is.na(x_num)        ~ 1L,
    x_num >= threshold  ~ 1L,
    TRUE                ~ 0L
  )
}

# Education: 5-level recode matching Stata (1/3=1, 4=2, 5=3, 6/7=4, 8/10=5)
code_education_5 <- function(x) {
  x <- to_na(x)
  dplyr::case_when(
    is.na(x)    ~ NA_integer_,
    x %in% 1:3  ~ 1L,
    x == 4      ~ 2L,
    x == 5      ~ 3L,
    x %in% 6:7  ~ 4L,
    x %in% 8:10 ~ 5L,
    TRUE        ~ NA_integer_
  )
}

code_edu_binary <- function(x, threshold = EDU_THRESHOLD_BIN) {
  x5 <- code_education_5(x)
  as.integer(!is.na(x5) & x5 >= threshold)
}

# Subjective wellbeing: s01 > 4
code_swb <- function(x, threshold = SWB_THRESHOLD) {
  x <- to_na(x)
  as.integer(!is.na(x) & x > threshold)
}

# Couple status: m36 in 1:3
code_couple <- function(x, codes = COUPLE_CODES) {
  x <- to_na(x)
  as.integer(!is.na(x) & x %in% codes)
}

# Institutional trust: keep as raw numeric
code_tinst <- function(x) to_na(x)

# Domain diversity: proportion of domains with membership
compute_diversity <- function(mat, cmax = N_DOMAINS) {
  rowSums(mat, na.rm = TRUE) / cmax
}

# Oster (2019) coefficient stability bound
oster_delta <- function(beta_r, beta_f, r2_r, r2_f, r2_max = NULL) {
  if (is.null(r2_max)) r2_max <- min(1.3 * r2_f, 1.0)
  if (abs(beta_r - beta_f) < 1e-10) return(Inf)
  (beta_f * (r2_max - r2_f)) / ((beta_r - beta_f) * (r2_f - r2_r))
}

ggsave_safe <- function(filename, plot, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(here::here("output", filename), plot,
                  width = width, height = height, dpi = dpi)
}

# ── Domain variables ──────────────────────────────────────────────────────────
# Order matches ISA Rmd 2021
DOMAIN_VARS <- c("nhg", "religious", "political", "union",
                 "professional", "charity", "sport", "student")

DOMAIN_LABELS <- c(
  nhg          = "Neighbors (JJVV)",
  religious    = "Religious",
  political    = "Political parties",
  union        = "Labor unions",
  professional = "Professional",
  charity      = "Charitable",
  sport        = "Sports clubs",
  student      = "Student org."
)

# ── Canonical STATE_MAP alpha / beta / gamma ──────────────────────────────────
# alpha = Isolate/Apathetic : low membership across all domains          clase=3 / ib3
# beta  = Closed/Clustering : concentrated memberships (mono-domain)     clase=1 / ib1
# gamma = Bridging/Broker   : diverse memberships (multi-domain)         clase=2 / ib2 [reference]
STATE_MAP <- tibble::tibble(
  greek       = c("\u03b1",      "\u03b2",      "\u03b3"),
  position    = c("alpha",       "beta",        "gamma"),
  clase       = c(3L,            1L,            2L),
  label_short = c("Isolate",     "Closed",      "Bridging"),
  label_long  = c(
    "\u03b1 (Isolate/Apathetic)",
    "\u03b2 (Closed/Clustering)",
    "\u03b3 (Bridging/Broker)"
  ),
  stata_ref   = c("ib3",         "ib1",         "ib2"),  # ib2 = reference in estimation.do
  is_ref      = c(FALSE,         FALSE,          TRUE)
)

POSITION_TO_CLASE <- setNames(STATE_MAP$clase, STATE_MAP$position)
CLASE_TO_LABEL    <- setNames(STATE_MAP$label_long, as.character(STATE_MAP$clase))

# ── Color palette — ggsci NPG (Nature Publishing Group) ────────
# Full NPG hex: #E64B35 #4DBBD5 #00A087 #3C5488 #F39B7F #8491B4 #91D1C2
STATE_PALETTE <- c(
  "\u03b1 (Isolate/Apathetic)"  = "#3C5488",   # dark navy   — withdrawn / low engagement
  "\u03b2 (Closed/Clustering)"  = "#4DBBD5",   # bright teal — focused cluster membership
  "\u03b3 (Bridging/Broker)"    = "#E64B35"    # warm red    — diverse cross-domain bridging
)

# Alias for scripts that use the shorter Greek label form
STATE_PALETTE_LEGACY <- c(
  "α (isolation)"  = "#3C5488",
  "β (clustering)" = "#4DBBD5",
  "γ (bridging)"   = "#E64B35"
)

# Wave / time-point colors (NPG, light → dark = early → late)
WAVE_PALETTE <- c(
  "2016" = "#91D1C2",   # mint  (lightest)
  "2018" = "#4DBBD5",   # teal
  "2022" = "#3C5488"    # dark navy (darkest)
)

STATE_LEVELS <- names(STATE_PALETTE)

# Labels for figures (keyed by numeric clase as character)
STATE_LABELS_FIG <- setNames(STATE_MAP$label_long, as.character(STATE_MAP$clase))
# Result: c("3"="\u03b1...", "1"="\u03b2...", "2"="\u03b3...")

# ── Plot themes (unified minimal style) ───────────────────────
theme_ssr <- function(base_size = 12, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Grid: horizontal only, very subtle
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.35),
      # Axes
      axis.title         = ggplot2::element_text(size = base_size * 0.92, color = "grey25"),
      axis.text          = ggplot2::element_text(size = base_size * 0.88, color = "grey30"),
      axis.ticks         = ggplot2::element_line(color = "grey80", linewidth = 0.3),
      axis.ticks.length  = ggplot2::unit(2.5, "pt"),
      # Facet strips — plain text, no background
      strip.text         = ggplot2::element_text(size = base_size * 0.95, color = "grey20",
                                                  margin = ggplot2::margin(b = 4)),
      strip.background   = ggplot2::element_blank(),
      # Legend
      legend.title       = ggplot2::element_blank(),
      legend.text        = ggplot2::element_text(size = base_size * 0.88, color = "grey20"),
      legend.position    = "top",
      legend.key.size    = ggplot2::unit(0.85, "lines"),
      legend.spacing.x   = ggplot2::unit(4, "pt"),
      # Panel border — consistent across all figures
      panel.border       = ggplot2::element_rect(color = "black", fill = NA,
                                                  linewidth = 0.65),
      # No titles or captions in any figure
      plot.title         = ggplot2::element_blank(),
      plot.subtitle      = ggplot2::element_blank(),
      plot.caption       = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(8, 12, 6, 8)
    )
}

theme_ssr_big <- function(base_size = 14, base_family = "sans") {
  theme_ssr(base_size = base_size, base_family = base_family)
}

# ── Startup summary ────────────────────────────────────────────────────────────────────────
message("[00_setup.R] loaded.")
message("  MEMBER_CODE_LOGIC = '", MEMBER_CODE_LOGIC, "' (any_member: c12 >= 2)")
message("  STATE_MAP: α=clase3/ib3 | β=clase1/ib1 | γ=clase2/ib2 [reference]")