# ==============================================================================
# any_member / 01_prepare.R
# Prepara el panel balanceado con codificación ANY_MEMBER (c12 >= 2)
# Output: output/dt_panel.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(here); library(panelr)
})

ROOT    <- here::here()
DIR_RAW <- file.path(ROOT, "data")
DIR_OUT <- file.path(ROOT, "any_member", "output")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# ── ÚNICA DIFERENCIA ENTRE PIPELINES ─────────────────────────────────────────
CODING <- "any_member"   # any_member/01_prepare.R  → c12 >= 2
# CODING <- "active_only" # active_only/01_prepare.R → c12 == 3
# ─────────────────────────────────────────────────────────────────────────────

ITEMS    <- c(paste0("c12_0", 1:8), "c12_09")
DOMAINS  <- c("nhg","religious","political","union",
               "professional","charity","sport","student","other")
WAVES    <- c(1L, 3L, 6L)
NA_CODES <- c(-999L, -888L, -777L, -666L)

cat("Pipeline:", CODING, "\n")
cat("Regla:   ", if (CODING=="any_member") "c12 >= 2 (inactivo o activo = 1)"
                 else                      "c12 == 3 (solo activo = 1)", "\n\n")

# ------------------------------------------------------------------------------
# 1. Cargar ELSOC
# ------------------------------------------------------------------------------
stopifnot(file.exists(file.path(DIR_RAW, "ELSOC_Long.RData")))
env   <- new.env()
nms   <- load(file.path(DIR_RAW, "ELSOC_Long.RData"), envir = env)
elsoc <- get(if ("elsoc_long_2016_2022" %in% nms) "elsoc_long_2016_2022" else nms[1], envir = env)

cat("ELSOC raw:", nrow(elsoc), "obs\n")

# ------------------------------------------------------------------------------
# 2. Filtros y selección de olas
# ------------------------------------------------------------------------------
if ("muestra"       %in% names(elsoc)) elsoc <- filter(elsoc, muestra == 1)
if ("tipo_atricion" %in% names(elsoc)) elsoc <- filter(elsoc, tipo_atricion == 1)

dt <- elsoc %>%
  filter(ola %in% WAVES) %>%
  mutate(ola = case_when(ola == 1L ~ 1L,
                         ola == 3L ~ 2L,
                         ola == 6L ~ 3L)) %>%
  mutate(across(where(is.numeric), ~ replace(.x, .x %in% NA_CODES, NA)))

cat("Después de filtros:", nrow(dt), "obs |", n_distinct(dt$idencuesta), "personas\n")

# ------------------------------------------------------------------------------
# 3. Codificación binaria de membresía
# ------------------------------------------------------------------------------
member_bin <- if (CODING == "any_member") {
  function(x) as.integer(!is.na(x) & x >= 2L)
} else {
  function(x) as.integer(!is.na(x) & x == 3L)
}

dt <- dt %>% mutate(across(all_of(ITEMS), member_bin))

# Conteo de dominios por persona-ola
dt <- dt %>%
  mutate(n_domains = rowSums(across(all_of(ITEMS)), na.rm = TRUE))

cat("\nDistribución n_domains (", CODING, "):\n")
print(table(dt$n_domains))

# ------------------------------------------------------------------------------
# 4. Variables de outcome y controles (según estimation.do)
# ------------------------------------------------------------------------------
dt <- dt %>%
  mutate(
    # Trust generalizada: c02 (1=sí, 2=depende→sí, 3=no)
    trust    = as.integer(c02 %in% c(1L, 3L)),     # 1=confía, 0=no confía

    # Trust vecinal: t01 > 3
    trust_nh = as.integer(!is.na(t01) & t01 > 3L),

    # Sociodemográficos
    woman    = as.integer(m0_sexo == 2L),
    age      = m0_edad,

    # Educación: 5 niveles → dummy high (técnica/university)
    education = case_when(
      m01 %in% 1:3 ~ 1L,   # básica
      m01 %in% 4:5 ~ 2L,   # media
      m01 %in% 6:7 ~ 3L,   # técnica
      m01 %in% 8:10~ 4L,   # universitaria
      TRUE          ~ NA_integer_
    ),
    edu_high = as.integer(!is.na(education) & education >= 3L),

    # Empleo: m02 <= 2
    employed = as.integer(!is.na(m02) & m02 <= 2L),

    # Bienestar subjetivo: s01 > 4
    swb      = as.integer(!is.na(s01) & s01 > 4L),

    # Pareja: m36 < 4
    couple   = as.integer(!is.na(m36) & m36 < 4L),

    # Confianza institucional (promedio de indicadores clave si existen)
    tinst    = case_when(
      !is.na(t02_01) ~ as.numeric(t02_01),
      TRUE            ~ NA_real_
    ),

    ola_f = factor(ola)
  )

# ------------------------------------------------------------------------------
# 5. Panel balanceado (3 olas completas en c12 items)
# ------------------------------------------------------------------------------
# Anclar a personas con 3 olas y sin NAs en items de membresía
ids_bal <- dt %>%
  group_by(idencuesta) %>%
  summarise(
    n_olas     = n_distinct(ola),
    any_na_items = any(is.na(rowSums(across(all_of(ITEMS))))),
    .groups = "drop"
  ) %>%
  filter(n_olas == 3L, !any_na_items) %>%
  pull(idencuesta)

dt_bal <- dt %>%
  filter(idencuesta %in% ids_bal) %>%
  select(idencuesta, ola, ola_f,
         all_of(ITEMS), n_domains,
         trust, trust_nh,
         woman, age, education, edu_high, employed, swb, couple, tinst)

cat("\nPanel balanceado:", n_distinct(dt_bal$idencuesta), "personas |",
    nrow(dt_bal), "obs\n")

cat("\nDistribución n_domains en panel balanceado:\n")
print(table(dt_bal$n_domains))

cat("\nResumen trust:\n")
cat("  trust    (gen):", round(mean(dt_bal$trust,    na.rm=TRUE)*100, 1), "%\n")
cat("  trust_nh (vec):", round(mean(dt_bal$trust_nh, na.rm=TRUE)*100, 1), "%\n")

# ------------------------------------------------------------------------------
# 6. Guardar
# ------------------------------------------------------------------------------
saveRDS(dt_bal, file.path(DIR_OUT, "dt_panel.rds"))
cat("\n✓ Guardado:", file.path(DIR_OUT, "dt_panel.rds"), "\n")
cat("  N =", n_distinct(dt_bal$idencuesta), "individuos ×",
    n_distinct(dt_bal$ola), "olas =", nrow(dt_bal), "obs\n")
cat("\n[01_prepare.R DONE]\n")
