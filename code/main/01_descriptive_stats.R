# ==============================================================================
# 01_descriptive_stats.R — Descriptivos y preparación del panel balanceado
# ------------------------------------------------------------------------------
# ACTUALIZADO post-comparación con estimation.do (Feb 2024).
# CAMBIOS:
#   (1) Filtros: muestra==1, tipo_atricion==1 (antes del filtro de olas)
#   (2) Trust vecinal: t01 (no c03)
#   (3) Membership: usa MEMBER_CODE_LOGIC de 00_setup.R
#   (4) Nuevas variables: swb, couple, t02_01 (controles del paper)
#   (5) Educación: recode 5 niveles igual que Stata
#   (6) Balanceo: replica lógica Stata (probit e(sample) + n_group==3)
#       → paso (a) drop obs con covariables missing;
#         paso (b) keep solo ids con 3 olas
#
# IMPORTANTE (estimation.do línea 76-82):
#   Stata balancea así:
#     probit social_trust ib2.clase i.m0_sexo m0_edad i.swb i.employed
#            i.couple i.t02_01 i.ola, vce(cluster idencuesta)
#     gen in_model = e(sample)
#     keep if in_model == 1
#   Esto elimina obs que le faltan CUALQUIER covariable del modelo completo.
#   Luego keep n_group == 3.
#   En R replicamos esto con complete.cases() sobre las mismas variables.
#
# Outputs:
#   output/summary_stats_balanced.csv
#   output/summary_stats_full.csv
#   output/membership_by_domain.csv
#   data/dt_analysis.rds   (dataset listo para LMM y modelos)
# ==============================================================================

source(here::here("code", "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(data.table)
})

# ==============================================================================
# 1. Cargar datos
# ==============================================================================
stop_if_missing(c(here::here("data", "ELSOC_Long.RData")))
obj_names <- load(here::here("data", "ELSOC_Long.RData"))

raw <- if ("elsoc_long_2016_2022" %in% obj_names) {
  get("elsoc_long_2016_2022")
} else if (length(obj_names) == 1) {
  get(obj_names[1])
} else {
  stop("No se puede identificar el objeto ELSOC. Disponibles: ",
       paste(obj_names, collapse = ", "))
}
dt_raw <- dplyr::as_tibble(raw)
message("Objeto cargado: ", nrow(dt_raw), " obs x ", ncol(dt_raw), " variables")

# ==============================================================================
# 2. FILTROS DE MUESTRA (replicando estimation.do líneas 65-67)
#    ANTES de filtrar por ola, para no distorsionar la composición
# ==============================================================================
dt <- dt_raw

if (FILTER_MUESTRA && "muestra" %in% names(dt)) {
  n_antes <- nrow(dt)
  dt <- dt %>% dplyr::filter(muestra == 1)
  message("Filtro muestra==1: ", n_antes, " → ", nrow(dt), " obs")
} else if (FILTER_MUESTRA) {
  warning("Variable 'muestra' no encontrada. Filtro no aplicado.")
}

if (FILTER_TIPO_ATRICION && "tipo_atricion" %in% names(dt)) {
  n_antes <- nrow(dt)
  dt <- dt %>% dplyr::filter(tipo_atricion == 1)
  message("Filtro tipo_atricion==1: ", n_antes, " → ", nrow(dt), " obs")
} else if (FILTER_TIPO_ATRICION) {
  warning("Variable 'tipo_atricion' no encontrada. Filtro no aplicado.")
}

# Filtrar olas de análisis
dt <- dt %>%
  dplyr::filter(ola %in% WAVES_RAW) %>%
  dplyr::mutate(
    id  = idencuesta,
    ola_rec = dplyr::case_when(
      ola == WAVES_RAW[1] ~ 1L,
      ola == WAVES_RAW[2] ~ 2L,
      ola == WAVES_RAW[3] ~ 3L
    ),
    wave_year = dplyr::case_when(
      ola_rec == 1 ~ WAVES_YEARS[1],
      ola_rec == 2 ~ WAVES_YEARS[2],
      ola_rec == 3 ~ WAVES_YEARS[3]
    )
  ) %>%
  dplyr::rename(ola_orig = ola, ola = ola_rec)

message("Obs en olas de análisis (1,3,6): ", nrow(dt))

# ==============================================================================
# 3. Membresías — usa MEMBER_CODE_LOGIC global
# ==============================================================================
missing_member <- setdiff(MEMBER_ITEMS, names(dt))
if (length(missing_member) > 0)
  stop("Ítems de membresía no encontrados: ", paste(missing_member, collapse=", "))

dt <- dt %>%
  dplyr::mutate(
    nhg          = member_binary(c12_01),
    religious    = member_binary(c12_02),
    sport        = member_binary(c12_03),   # c12_03 = deportiva
    charity      = member_binary(c12_04),   # c12_04 = caridad
    political    = member_binary(c12_05),   # c12_05 = partido/mov
    professional = member_binary(c12_06),
    union        = member_binary(c12_07),
    student      = member_binary(c12_08)
  )

# Verificar si existe c12_09 (ítem 9, "otra", usado en ISA Rmd 2021)
if (MEMBER_ITEMS_9TH %in% names(dt)) {
  dt <- dt %>% dplyr::mutate(otra = member_binary(c12_09))
  if (!"otra" %in% DOMAIN_VARS) DOMAIN_VARS <<- c(DOMAIN_VARS, "otra")
  N_DOMAINS <<- 9L
  message("c12_09 encontrado y agregado como 9vo dominio.")
} else {
  message("c12_09 no encontrado. Usando 8 dominios.")
}

# ==============================================================================
# 4. Variables de outcome
# ==============================================================================
for (v in c(TRUST_GEN_VAR, TRUST_NH_VAR)) {
  if (!v %in% names(dt)) stop("Variable no encontrada: ", v,
                               "\n¿Es t01 para trust vecinal? Verificar codebook.")
}

dt <- dt %>%
  dplyr::mutate(
    trust    = code_gen_trust(.data[[TRUST_GEN_VAR]]),
    trust_nh = code_nh_trust(.data[[TRUST_NH_VAR]])
  )

message("Trust generalizada: ", round(mean(dt$trust, na.rm=TRUE), 3),
        " | Trust vecinal: ", round(mean(dt$trust_nh, na.rm=TRUE), 3))

# ==============================================================================
# 5. Covariables (replicando estimation.do)
# ==============================================================================
dt <- dt %>%
  dplyr::mutate(
    # Sexo: m0_sexo (1=hombre, 2=mujer en ELSOC) → woman = 0/1
    woman     = dplyr::case_when(m0_sexo == 1 ~ 0L, m0_sexo == 2 ~ 1L, TRUE ~ NA_integer_),
    edad      = to_na(m0_edad),
    education = code_education_5(.data[[EDU_VAR]]),  # 5 niveles como Stata
    edu_bin   = code_edu_binary(.data[[EDU_VAR]]),

    # NUEVA: employed = 1 if m02 <= 2 (replicando estimation.do exacto)
    employed  = if (EMPLOY_VAR %in% names(.)) {
      code_employment(.data[[EMPLOY_VAR]])
    } else {
      warning("'", EMPLOY_VAR, "' no encontrado."); NA_integer_
    },

    # NUEVA: swb = 1 if s01 > 4 (estimation.do línea 62)
    swb = if (SWB_VAR %in% names(.)) {
      code_swb(.data[[SWB_VAR]])
    } else {
      warning("'", SWB_VAR, "' no encontrado."); NA_integer_
    },

    # NUEVA: couple = 1 if m36 < 4 (estimation.do línea 52)
    couple = if (COUPLE_VAR %in% names(.)) {
      code_couple(.data[[COUPLE_VAR]])
    } else {
      warning("'", COUPLE_VAR, "' no encontrado."); NA_integer_
    },

    # NUEVA: t02_01 (confianza institucional — control en modelos Stata)
    tinst = if (TINST_VAR %in% names(.)) {
      code_tinst(.data[[TINST_VAR]])
    } else {
      warning("'", TINST_VAR, "' no encontrado."); NA_real_
    }
  )

# ==============================================================================
# 6. Variables derivadas de membresía
# ==============================================================================
dt <- dt %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    membership_count = sum(dplyr::c_across(dplyr::all_of(DOMAIN_VARS)), na.rm = TRUE),
    domain_diversity = membership_count / N_DOMAINS,
    multi_domain     = as.integer(membership_count >= 2L)
  ) %>%
  dplyr::ungroup()

# ==============================================================================
# 7. BALANCEO DEL PANEL — replicando lógica Stata
#
#   Stata (estimation.do líneas 76-87):
#     (a) Probit con todos los controles → gen in_model = e(sample)
#         → equivale a: obs con NO-missing en TODAS las variables del modelo
#     (b) keep n_group == 3 (solo ids con 3 olas)
# ==============================================================================

# Paso (a): variables que Stata usa para definir la muestra analítica
stata_model_vars <- c("trust","swb","employed","couple","tinst","woman","edad")
stata_model_vars <- intersect(stata_model_vars, names(dt))
stata_model_vars <- stata_model_vars[sapply(stata_model_vars,
                    function(v) any(!is.na(dt[[v]])))]

# Obs con complete data en todas las variables del modelo
complete_mask <- complete.cases(dt[, stata_model_vars])
dt$in_model   <- complete_mask

message("Paso (a) — complete cases en vars. del modelo: ",
        sum(complete_mask), " de ", nrow(dt), " obs (",
        round(mean(complete_mask)*100, 1), "%)")

# Paso (b): ids con exactamente 3 olas (después de filtrar complete cases)
panel_counts <- dt %>%
  dplyr::filter(in_model) %>%
  dplyr::count(id)

balanced_ids <- panel_counts %>% dplyr::filter(n == 3) %>% dplyr::pull(id)

dt$in_balanced <- dt$id %in% balanced_ids & dt$in_model

n_balanced <- length(balanced_ids)
n_total    <- length(unique(dt$id))
message("Paso (b) — Panel balanceado: n = ", n_balanced, " individuos (",
        round(n_balanced / n_total * 100, 1), "% de retención)")

# ==============================================================================
# 8. Estadísticas descriptivas
# ==============================================================================
summary_fun <- function(d) {
  d %>%
    dplyr::group_by(ola) %>%
    dplyr::summarise(
      N                  = dplyr::n(),
      age_mean           = mean(edad,             na.rm = TRUE),
      age_sd             = sd(edad,              na.rm = TRUE),
      pct_woman          = mean(woman,            na.rm = TRUE),
      pct_edu_secondary  = mean(edu_bin,          na.rm = TRUE),
      pct_employed       = mean(employed,         na.rm = TRUE),
      pct_swb            = mean(swb,              na.rm = TRUE),
      pct_couple         = mean(couple,           na.rm = TRUE),
      generalized_trust  = mean(trust,            na.rm = TRUE),
      neighborhood_trust = mean(trust_nh,         na.rm = TRUE),
      membership_mean    = mean(membership_count, na.rm = TRUE),
      membership_sd      = sd(membership_count,  na.rm = TRUE),
      pct_multi_domain   = mean(multi_domain,     na.rm = TRUE),
      domain_diversity   = mean(domain_diversity, na.rm = TRUE),
      .groups = "drop"
    )
}

stats_balanced <- dt %>% dplyr::filter(in_balanced) %>% summary_fun()
stats_full     <- dt %>% summary_fun()

readr::write_csv(stats_balanced, here::here("output", "summary_stats_balanced.csv"))
readr::write_csv(stats_full,     here::here("output", "summary_stats_full.csv"))

cat("\n--- Estadísticas (panel balanceado) ---\n")
print(stats_balanced)

# ==============================================================================
# 9. Membresía por dominio y ola
# ==============================================================================
domain_tbl <- dt %>%
  dplyr::filter(in_balanced) %>%
  dplyr::group_by(ola) %>%
  dplyr::summarise(
    N = dplyr::n(),
    dplyr::across(dplyr::all_of(DOMAIN_VARS), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(dplyr::all_of(DOMAIN_VARS),
                      names_to = "domain", values_to = "pct_member") %>%
  dplyr::arrange(ola, dplyr::desc(pct_member))

readr::write_csv(domain_tbl, here::here("output", "membership_by_domain.csv"))

# ==============================================================================
# 10. Verificación de distribuciones (diagnósticos)
# ==============================================================================
cat("\n--- DIAGNÓSTICO: Trust generalizada (c02) en panel balanceado, ola 1 ---\n")
dt %>% dplyr::filter(in_balanced, ola == 1) %>%
  dplyr::count(.data[[TRUST_GEN_VAR]], trust) %>% print()

cat("\n--- DIAGNÓSTICO: Trust vecinal (", TRUST_NH_VAR, ") en panel balanceado, ola 1 ---\n")
dt %>% dplyr::filter(in_balanced, ola == 1) %>%
  dplyr::count(.data[[TRUST_NH_VAR]], trust_nh) %>% print()

cat("\n--- DIAGNÓSTICO: swb en panel balanceado, ola 1 ---\n")
dt %>% dplyr::filter(in_balanced, ola == 1) %>%
  dplyr::count(.data[[SWB_VAR]], swb) %>% print()

cat("\n--- DIAGNÓSTICO: Membresía (", MEMBER_CODE_LOGIC, ") ola 1 ---\n")
cat("Proporción con >= 1 membresía:",
    round(mean(dt$membership_count[dt$in_balanced & dt$ola==1] > 0, na.rm=TRUE), 3), "\n")

# ==============================================================================
# 11. Guardar dataset listo para análisis
# ==============================================================================
dt_analysis <- dt %>%
  dplyr::filter(in_balanced) %>%
  dplyr::select(
    id, ola, wave_year, trust, trust_nh, woman, edad, education, edu_bin,
    employed, swb, couple, tinst,
    dplyr::all_of(DOMAIN_VARS),
    membership_count, domain_diversity, multi_domain
  )

saveRDS(dt_analysis, here::here("data", "dt_analysis.rds"))
message("\n[01_descriptive_stats.R] Listo.")
message("  Panel balanceado: n = ", length(balanced_ids),
        " individuos | ", nrow(dt_analysis), " obs person-wave")
message("  Saved: data/dt_analysis.rds")
message("\n  RECORDATORIO: Si trust_nh da 0% o 100%, verificar que t01 existe en el dataset.")
message("  RECORDATORIO: Si employed es todo NA, verificar que m02 existe.")
message("  RECORDATORIO: MEMBER_CODE_LOGIC = '", MEMBER_CODE_LOGIC, "'")
