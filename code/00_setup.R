# ==============================================================================
# 00_setup.R — Global setup (packages, paths, helpers)
# ------------------------------------------------------------------------------
# VERSIÓN ACTUALIZADA post-comparación con:
#   estimation.do (Feb 2024), estimation_2.do (Feb 2024)
#   membershipspatterns_trustISA2021.Rmd (ISA Forum 2021)
#
# CAMBIOS CRÍTICOS vs. versión anterior:
#   (1) TRUST_NH_VAR: "c03" → "t01"       [estimation.do líneas 35-36]
#   (2) to_na(): agrega -666               [estimation.do mvdecode -666]
#   (3) MEMBER_CODE_LOGIC: flag explícito  [ISA Rmd usa c12 >= 2]
#   (4) Filtros muestra==1, tipo_atricion==1 [estimation.do líneas 65-66]
#   (5) Controles nuevos: swb, couple, t02_01 [estimation.do líneas 62-77]
#   (6) EDU recode de 5 niveles            [estimation.do línea 41]
#   (7) STATE_MAP canónico α/β/γ           [coherencia paper-código]
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(data.table)
})

DIR_DATA <- here::here("data")
DIR_OUT  <- here::here("output")
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# ==============================================================================
# GLOBAL KNOBS
# ==============================================================================

WAVES_RAW   <- c(1, 3, 6)
WAVES_YEARS <- c(2016, 2018, 2022)

MEMBER_ITEMS <- paste0("c12_0", 1:8)
# c12_09 ("otra organización") se incluía en ISA Rmd 2021.
# Verificar si existe en ELSOC Long 2016-2022 y agregar si corresponde.
MEMBER_ITEMS_9TH <- "c12_09"   # poner NULL para excluir

# ==============================================================================
# MEMBERSHIP CODING — VERIFICAR CONTRA CODEBOOK ELSOC
# ------------------------------------------------------------------------------
# ELSOC Long 2016-2022 — CODEBOOK CONFIRMADO DEL DATASET REAL:
#   1 = No es miembro
#   2 = Miembro inactivo
#   3 = Miembro activo
#
# any_member  = c12 >= 2  (inactivo o activo — usado en ISA Rmd 2021)
# active_only = c12 == 3  (solo activo — definición estricta)
# ==============================================================================
MEMBER_CODE_LOGIC <- "any_member"
# "any_member"  → c12 >= 2  (miembro inactivo o activo) ← CORRECTO según ELSOC [ANÁLISIS PRINCIPAL]
# "active_only" → c12 == 3  (solo miembro activo, definición estricta)          [SENSIBILIDAD S7]

# TRUST — generalizada (confirmado: estimation.do línea 31)
TRUST_GEN_VAR  <- "c02"
TRUST_GEN_CODE <- 1L

# TRUST — vecinal
# FIX: estimation.do usa t01, NO c03 (error crítico en versión previa)
TRUST_NH_VAR       <- "t01"
TRUST_NH_THRESHOLD <- 4L    # t01 > 3 ↔ t01 >= 4

# EMPLEO — estimation.do líneas 47-49: employed = 1 if m02 <= 2
EMPLOY_VAR   <- "m02"
EMPLOY_CODES <- 1L:2L

# EDUCACIÓN — recode 5 niveles (estimation.do línea 41)
EDU_VAR <- "m01"
EDU_THRESHOLD_BIN <- 4L   # >= 4 en escala de 5 = media completa o más

# BIENESTAR SUBJETIVO — estimation.do líneas 62-63: swb = 1 if s01 > 4
SWB_VAR       <- "s01"
SWB_THRESHOLD <- 4L

# PAREJA — estimation.do líneas 52-53: couple = 1 if m36 < 4
COUPLE_VAR   <- "m36"
COUPLE_CODES <- 1L:3L

# CONFIANZA INSTITUCIONAL (control en modelos Stata)
TINST_VAR <- "t02_01"

# FILTROS DE MUESTRA — estimation.do líneas 65-66
FILTER_MUESTRA       <- TRUE   # keep if muestra == 1
FILTER_TIPO_ATRICION <- TRUE   # keep if tipo_atricion == 1

# ==============================================================================
# CONTROL VARIABLES — definidos una sola vez, heredados por todos los scripts
# Replicando: estimation.do L.101-102
#   xtprobit social_trust ib2.clase i.m0_sexo m0_edad i.swb
#             i.employed i.couple i.t02_01 i.ola
# ==============================================================================
CONTROLS_BASE  <- c("edad", "woman", "education")
CONTROLS_EXTRA <- c("swb", "employed", "couple", "tinst")
CONTROLS_FULL  <- c(CONTROLS_BASE, CONTROLS_EXTRA)

# LMM
K_GRID     <- 1:5
K_BASELINE <- 3
N_STARTS   <- 30
MOD_TRANS  <- 1
TOL        <- 1e-8
MAXIT      <- 1000

N_DOMAINS <- 8L

# ==============================================================================
# HELPERS
# ==============================================================================

stop_if_missing <- function(paths) {
  miss <- paths[!file.exists(paths)]
  if (length(miss) > 0)
    stop("Missing files:\n- ", paste(miss, collapse = "\n- "))
}

# ELSOC missing codes: -999=No responde, -888=No sabe, -777=Error técnico, -666=Encuesta incompleta
to_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% c(-999, -888, -777, -666, 777777, 999999)] <- NA_real_
  x
}

# Membership: dos lógicas según codebook ELSOC
member_binary <- function(x, logic = MEMBER_CODE_LOGIC) {
  x <- to_na(x)
  if (logic == "any_member") {
    as.integer(!is.na(x) & x >= 2L)
  } else {
    as.integer(!is.na(x) & x == 3L)
  }
}
member_active <- function(x) member_binary(x)

# Employment: estimation.do m02 <= 2
code_employment <- function(x, codes = EMPLOY_CODES) {
  x <- to_na(x)
  as.integer(!is.na(x) & x %in% codes)
}

# Trust generalizada — replica EXACTAMENTE estimation.do:
#   mvdecode _all, mv(-888/-999/-666)     → missing en Stata
#   gen social_trust = 1 if c02 == 1
#   replace social_trust = 0 if c02 > 1  ← EN STATA: missing(.) > 1 es TRUE
#
# Por eso individuos con c02==-888/-999 reciben social_trust=0 en Stata
# (no missing). En R debemos replicar esto NO llamando to_na() primero,
# sino mapeando los valores de missing codes directamente a 0.
code_gen_trust <- function(x, trust_code = TRUST_GEN_CODE) {
  x_num <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    x_num == trust_code                          ~ 1L,  # c02==1 -> trust
    !is.na(x_num) & x_num != trust_code          ~ 0L,  # c02==2,3 -> no trust
    is.na(x_num)                                 ~ 0L   # c02==missing -> 0
    # (Stata: missing > 1 es TRUE, entonces replace social_trust=0 aplica)
  )
}

# Trust vecinal — replica estimation.do:
#   gen trust_nb = 1 if t01 > 3
#   replace trust_nb = 0 if t01 <= 3
# En Stata: missing(.) <= 3 es FALSE (. = +inf), entonces t01==missing -> trust_nb=1.
# Pero en la práctica t01 missing es raro y el efecto en N es pequeño.
# Por seguridad: missing -> 1 (igual que Stata).
code_nh_trust <- function(x, threshold = TRUST_NH_THRESHOLD) {
  x_num <- suppressWarnings(as.numeric(x))
  dplyr::case_when(
    is.na(x_num)        ~ 1L,  # missing -> 1 (Stata: . <= 3 es FALSE)
    x_num >= threshold  ~ 1L,  # t01 >= 4 -> trust
    TRUE                ~ 0L   # t01 <= 3 -> no trust
  )
}

# Educación 5 niveles (estimation.do: recode m01 1/3=1, 4=2, 5=3, 6/7=4, 8/10=5)
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

# SWB: s01 > 4 (estimation.do)
code_swb <- function(x, threshold = SWB_THRESHOLD) {
  x <- to_na(x)
  as.integer(!is.na(x) & x > threshold)
}

# Couple: m36 < 4 (estimation.do)
code_couple <- function(x, codes = COUPLE_CODES) {
  x <- to_na(x)
  as.integer(!is.na(x) & x %in% codes)
}

# Confianza institucional: mantener como raw
code_tinst <- function(x) to_na(x)

# Domain diversity
compute_diversity <- function(mat, cmax = N_DOMAINS) {
  rowSums(mat, na.rm = TRUE) / cmax
}

# Oster (2019)
oster_delta <- function(beta_r, beta_f, r2_r, r2_f, r2_max = NULL) {
  if (is.null(r2_max)) r2_max <- min(1.3 * r2_f, 1.0)
  if (abs(beta_r - beta_f) < 1e-10) return(Inf)
  (beta_f * (r2_max - r2_f)) / ((beta_r - beta_f) * (r2_f - r2_r))
}

ggsave_safe <- function(filename, plot, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(here::here("output", filename), plot,
                  width = width, height = height, dpi = dpi)
}

# ==============================================================================
# DOMAIN VARIABLES
# Orden igual que ISA Rmd 2021
# ==============================================================================
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

# ==============================================================================
# STATE MAP CANÓNICO α / β / γ
# ------------------------------------------------------------------------------
# Este es el mapeo definitivo entre:
#   (a) símbolo griego del paper,
#   (b) nombre sustantivo (paper text),
#   (c) label en código R (position en dt_states),
#   (d) clase numérica (para ib1/ib2 en Stata),
#   (e) referencia en estimation.do.
#
# REGLA: NO alterar este bloque sin actualizar paper, tablas y figuras.
#
# α  =  Apathetic / Isolate:  baja membresía en todos los dominios
# β  =  Closed / Clustering:  membresías concentradas (mono-dominio, alta densidad)
# γ  =  Bridging / Broker:    membresías diversas (multi-dominio) ← REFERENCIA
# ==============================================================================
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
  stata_ref   = c("ib3",         "ib1",         "ib2"),  # ib2 = ref en estimation.do
  is_ref      = c(FALSE,         FALSE,          TRUE)   # γ = categoría de referencia
)

# Lookup rápido: position → clase numérica
POSITION_TO_CLASE <- setNames(STATE_MAP$clase, STATE_MAP$position)
# Lookup rápido: clase numérica → label_long
CLASE_TO_LABEL    <- setNames(STATE_MAP$label_long, as.character(STATE_MAP$clase))

# ==============================================================================
# PALETA DE COLORES — consistente con STATE_MAP
# ==============================================================================
STATE_PALETTE <- c(
  "\u03b1 (Isolate/Apathetic)"  = "#59A14F",   # verde
  "\u03b2 (Closed/Clustering)"  = "#4E79A7",   # azul
  "\u03b3 (Bridging/Broker)"    = "#E15759"    # rojo ← referencia
)

# Alias para compatibilidad con código anterior (sin modificar 002)
STATE_PALETTE_LEGACY <- c(
  "α (isolation)"  = "#59A14F",
  "β (clustering)" = "#4E79A7",
  "γ (bridging)"   = "#E15759"
)

STATE_LEVELS <- names(STATE_PALETTE)

# Labels para figuras (keyed by clase numérica como character)
STATE_LABELS_FIG <- setNames(STATE_MAP$label_long, as.character(STATE_MAP$clase))
# Resultado: c("3"="\u03b1...", "1"="\u03b2...", "2"="\u03b3...")

# ==============================================================================
# THEMES
# ==============================================================================
theme_ssr <- function(base_size = 13, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title         = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold"),
      legend.title       = ggplot2::element_blank(),
      legend.position    = "top"
    )
}

theme_ssr_big <- function(base_size = 16, base_family = "sans") {
  theme_ssr(base_size = base_size, base_family = base_family) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size))
}

# ==============================================================================
# MENSAJES DE VERIFICACIÓN
# ==============================================================================
message("[00_setup.R] ACTUALIZADO — CAMBIOS CRÍTICOS:")
message("  (1) TRUST_NH_VAR = 't01' (era 'c03')")
message("  (2) to_na() incluye -999/-888/-777/-666")
message("  (3) MEMBER_CODE_LOGIC = '", MEMBER_CODE_LOGIC, "'")
message("  (4) Filtros: muestra==1, tipo_atricion==1")
message("  (5) Controles nuevos: swb (s01>4), couple (m36<4), t02_01")
message("  (7) STATE_MAP canónico definido:")
message("      \u03b1 (alpha) = Isolate/Apathetic | clase 3")
message("      \u03b2 (beta)  = Closed/Clustering  | clase 1")
message("      \u03b3 (gamma) = Bridging/Broker    | clase 2 [REFERENCIA ib2]")
message("")
message("  CODEBOOK CONFIRMADO (ELSOC Long 2016-2022):")
message("  c12: 1=No es miembro | 2=Miembro inactivo | 3=Miembro activo")
message("  MEMBER_CODE_LOGIC='any_member' (c12>=2) = inactivo o activo \u2713")