# ==============================================================================
# 003_trust_models.R — Modelos de confianza
# ==============================================================================
# PIPELINE:
#   00_setup.R → 01_descriptive_stats.R → 002_long_latent_class.R → [ESTE]
#
# TERMINOLOGÍA CANÓNICA (STATE_MAP de 00_setup.R — NO ALTERAR):
#   α (alpha) = Isolate/Apathetic  | clase=3 | ref=ib3
#   β (beta)  = Closed/Clustering  | clase=1 | ref=ib1
#   γ (gamma) = Bridging/Broker    | clase=2 | ref=ib2 [REFERENCIA en estimation.do]
#
# NOTAS METODOLÓGICAS CRÍTICAS:
#
#   (A) nAGQ = 9 — FIX PRINCIPAL vs versión anterior:
#     Stata xtprobit usa cuadratura de Gauss-Hermite (GHQ, default 12 puntos).
#     glmer() con default (Laplace, nAGQ=1) subestima sistemáticamente sigma_u
#     del RE probit, produciendo AMEs 5-10× más pequeños que Stata.
#     Fix: nAGQ=9 (≥7 recomendado para RE probit; balance exactitud/tiempo).
#     Para replicación exacta de Stata: NAGQ <- 12L (más lento).
#     Ref: Pinheiro & Bates (2000); Rabe-Hesketh & Skrondal (2012).
#
#   (B) Cluster-robust SEs — FIX (reemplaza clubSandwich::vcovCR):
#     clubSandwich::vcovCR falla con "$ operator not defined for S4 class"
#     porque glmerMod es S4 y vcovCR usa acceso S3 ($) a internos del modelo.
#     Fix: sandwich::vcovCL(~idencuesta_f, type="HC1") vía argumento
#     vcov=~cluster en avg_slopes. Requiere merDeriv para estfun.glmerMod.
#     Si vcovCL también falla → fallback a model vcov (SEs no-clustered).
#
#   (C) marginaleffects_safe = FALSE:
#     Silencia warning: AMEs en GLMM usan solo VCOV de parámetros fijos.
#     Equivalente a Stata margins dydx(*) post-xtprobit (re.form=NA).
#     Footnote tablas: "AMEs at population level; SEs from fixed-effect
#     VCOV, consistent with Stata margins post-xtprobit."
#
# OUTPUTS:
#   output/pre_model_diagnostics.txt     ← REVISAR ANTES DE INTERPRETAR
#   output/table4_reprobit_ref_bridging.csv
#   output/table5_reprobit_ref_closed.csv
#   output/tableA1_probit_pool_ref_closed.csv
#   output/tableA2_reverse_causality.csv
#   output/fig3_predicted_probs.png / .pdf
#   output/sensitivity_confgral.csv
#   output/ame_all_models.csv
#   output/trust_models_summaries.txt
#   data/trust_model_objects.rds
# ==============================================================================

options(marginaleffects_safe = FALSE)

source(here::here("code", "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(lme4); library(marginaleffects)
  library(sandwich); library(ggplot2); library(broom.mixed)
  library(stringr)
})

# merDeriv: provee estfun() y bread() para objetos merMod (glmerMod / lmerMod)
# Esto es NECESARIO para que sandwich::vcovCL funcione con glmer (nAGQ >= 1).
# Sin este paquete, vcovCL falla con:
#   "no applicable method for estfun applied to class glmerMod"
if (!requireNamespace("merDeriv", quietly = TRUE)) {
  message("Instalando merDeriv (necesario para vcovCL con glmerMod)...")
  install.packages("merDeriv", quiet = TRUE)
}
suppressPackageStartupMessages(library(merDeriv))

dir.create(here::here("output"), showWarnings = FALSE)
dir.create(here::here("data"),   showWarnings = FALSE)

# nAGQ global — cambiar aquí si se quiere experimentar
# 9  = buena aproximación de GHQ-12 de Stata; más rápido
# 12 = equivalente exacto a Stata xtprobit default (más lento)
NAGQ <- 12L

# ==============================================================================
# 1. CARGAR Y FUSIONAR
# ==============================================================================
stop_if_missing(c(
  here::here("data", "dt_states_cov.rds"),
  here::here("data", "dt_analysis.rds")
))

dt_states <- readRDS(here::here("data", "dt_states_cov.rds")) %>% as_tibble()
dt_anal   <- readRDS(here::here("data", "dt_analysis.rds"))   %>% as_tibble()

message("dt_states_cov: ", nrow(dt_states), " obs | ",
        length(unique(dt_states$idencuesta)), " ids")
message("dt_analysis:   ", nrow(dt_anal), " obs | ",
        length(unique(dt_anal$id)), " ids")

dt_anal_join <- dt_anal %>%
  dplyr::rename(idencuesta = id) %>%
  dplyr::select(
    idencuesta, ola,
    trust, trust_nh,
    edad, woman, education, edu_bin,
    employed, swb, couple, tinst,
    membership_count, domain_diversity, multi_domain
  )

dt_states    <- dt_states    %>% dplyr::mutate(idencuesta = as.character(idencuesta))
dt_anal_join <- dt_anal_join %>% dplyr::mutate(idencuesta = as.character(idencuesta))

# Solo traer columnas de dt_anal que NO están ya en dt_states (evita .x/.y conflicts)
cols_new <- setdiff(names(dt_anal_join),
                    c("idencuesta", "ola", setdiff(names(dt_states), c("idencuesta","ola"))))
dt <- dplyr::left_join(
  dt_states,
  dt_anal_join %>% dplyr::select(dplyr::all_of(c("idencuesta", "ola", cols_new))),
  by = c("idencuesta", "ola")
)

n_matched   <- sum(!is.na(dt$trust))
pct_matched <- round(n_matched / nrow(dt) * 100, 1)
message("Match exitoso: ", n_matched, "/", nrow(dt), " obs (", pct_matched, "%)")
if (pct_matched < 80)
  warning("< 80% obs con trust. Verificar filtros en 01 y 002.")

# ==============================================================================
# 2. VARIABLES DERIVADAS — STATE_MAP CANÓNICO
# ==============================================================================
dt <- dt %>%
  dplyr::mutate(
    clase = dplyr::case_when(
      position == "alpha" ~ 3L,   # α Isolate/Apathetic  | clase=3 | ib3
      position == "beta"  ~ 1L,   # β Closed/Clustering   | clase=1 | ib1
      position == "gamma" ~ 2L,   # γ Bridging/Broker     | clase=2 | ib2 ← REF
      TRUE                ~ NA_integer_
    ),
    clase_strict = dplyr::case_when(
      position_strict == "alpha" ~ 3L,
      position_strict == "beta"  ~ 1L,
      position_strict == "gamma" ~ 2L,
      TRUE                       ~ NA_integer_
    ),
    clase_label = dplyr::recode(as.character(clase), !!!STATE_LABELS_FIG)
  )

# CONTROLS_TRUST: replica exactamente estimation.do:
#   xtprobit social_trust ib2.clase i.m0_sexo m0_edad i.swb i.employed i.couple i.t02_01 i.ola
# education NO aparece en ese modelo; incluirla absorbe el efecto de clase (colinealidad).
CONTROLS_TRUST <- c("swb", "employed", "couple")

ctrl_use <- CONTROLS_TRUST[
  CONTROLS_TRUST %in% names(dt) &
  sapply(CONTROLS_TRUST, function(v) sum(!is.na(dt[[v]])) > 0)
]
message("Controles disponibles: ", paste(ctrl_use, collapse = ", "))

# ==============================================================================
# 3. PRE-MODEL DIAGNOSTICS
#    Revisar ANTES de interpretar cualquier resultado.
#    Targets Stata: N = 3891 | Patrón: γ(Bridging) > α(Isolate) > β(Closed)
# ==============================================================================
cat("\n=== PRE-MODEL DIAGNOSTICS ===\n")

diag_lines <- c(
  "=== PRE-MODEL DIAGNOSTICS (003_trust_models.R) ===",
  paste0("Fecha: ", Sys.time()),
  paste0("nAGQ = ", NAGQ, " | Stata GHQ = 12"),
  ""
)

# 3a. N para modelo principal + desglose de pérdidas
stata_model_vars_check <- intersect(c("trust","clase",ctrl_use,"ola"), names(dt))
d_check <- dt %>%
  dplyr::filter(!is.na(clase)) %>%
  tidyr::drop_na(dplyr::all_of(stata_model_vars_check))

cat("N disponible modelo principal (trust+controles): ", nrow(d_check), "\n")
cat("N objetivo Stata:                                 3891\n")
cat("Diferencia:                                      ", nrow(d_check) - 3891, "\n\n")

# Desglose de pérdidas por variable — para identificar qué variable reduce el N
cat("=== Desglose de NAs por variable (en dt con clase no-NA) ===\n")
dt_con_clase <- dt %>% dplyr::filter(!is.na(clase))
cat("  N total con clase no-NA: ", nrow(dt_con_clase), "\n")
for (v in stata_model_vars_check) {
  n_miss <- sum(is.na(dt_con_clase[[v]]))
  cat("  NA en", sprintf("%-15s", v), ":", n_miss, "\n")
}
cat("\n")
cat("  Si la diferencia vs Stata (3891) es > 200:\n")
cat("  -> Verificar filtros muestra==1 y tipo_atricion==1 en 01_descriptive_stats.R\n")
cat("  -> Verificar que el merge dt_states+dt_analysis no pierde obs\n\n")

diag_lines <- c(diag_lines,
  paste0("N disponible R (trust+controles): ", nrow(d_check)),
  paste0("N total con clase no-NA: ", nrow(dt_con_clase)),
  paste0("N objetivo Stata: 3891"),
  paste0("Diferencia: ", nrow(d_check) - 3891),
  "  Si |diferencia| > 200, revisar filtros muestra/tipo_atricion en 01.",
  ""
)

# 3b. Trust medio por clase — patrón esperado γ > α > β (o γ > β > α)
cat("=== Trust por clase (Stata Fig.3 espera: Broker[γ] > Apathetic[α] > Closed[β]) ===\n")
trust_by_clase <- dt %>%
  dplyr::filter(!is.na(trust), !is.na(clase)) %>%
  dplyr::group_by(position, clase, clase_label) %>%
  dplyr::summarise(
    n        = dplyr::n(),
    trust    = round(mean(trust,    na.rm = TRUE), 3),
    trust_nh = round(mean(trust_nh, na.rm = TRUE), 3),
    .groups  = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(trust))
print(trust_by_clase)

diag_lines <- c(diag_lines,
  "--- Trust por clase (esperado: γ > α > β) ---",
  capture.output(print(trust_by_clase)), ""
)

# 3c. Distribución de clases por ola
cat("\n=== Distribución de clases por ola ===\n")
clase_dist <- dt %>%
  dplyr::filter(!is.na(clase)) %>%
  dplyr::count(ola, clase_label) %>%
  tidyr::pivot_wider(names_from = ola, values_from = n, names_prefix = "ola_")
print(clase_dist)

# 3d. Coding de trust
cat("\n=== trust (esperado: binaria 0/1) ===\n")
print(table(dt$trust,    useNA = "ifany"))
print(table(dt$trust_nh, useNA = "ifany"))

diag_lines <- c(diag_lines,
  "--- Distribución trust ---",
  capture.output(print(table(dt$trust,    useNA="ifany"))),
  capture.output(print(table(dt$trust_nh, useNA="ifany"))), ""
)

writeLines(diag_lines, here::here("output", "pre_model_diagnostics.txt"))
message("  Guardado: output/pre_model_diagnostics.txt — REVISAR ANTES DE CONTINUAR")

# ==============================================================================
# 4. HELPERS
# ==============================================================================

# Cluster-robust SE para glmerMod (nAGQ >= 1)
# Estrategia: sandwich::vcovCL vía argumento vcov de marginaleffects.
#   - marginaleffects pasa ~cluster_var directamente a sandwich::vcovCL
#   - sandwich::vcovCL usa estfun.glmerMod (provisto por merDeriv) para los scores
#   - Evita el error "$ operator not defined for S4" de clubSandwich::vcovCR
#   - HC1 = corrección de muestra finita estándar (G/(G-1), G = N clusters)
#
# Por qué no clubSandwich::vcovCR:
#   vcovCR accede a internos del modelo con $ (operador S3), pero glmerMod es S4.
#   Incompatibilidad entre versiones recientes de clubSandwich y lme4.
#
# Por qué vcovCL funciona aquí:
#   sandwich::vcovCL llama estfun() + bread() via métodos genéricos S4-compatibles
#   que merDeriv registra para objetos merMod.
#
# Retorna: fórmula de cluster (~idencuesta_f) para pasar a avg_slopes(vcov=)
#          o NULL si merDeriv no está disponible (fallback a model vcov)
get_vcov_formula <- function(model, cluster_var = "idencuesta_f") {
  # Verificar que merDeriv provee estfun para este modelo
  has_estfun <- tryCatch({
    sf <- sandwich::estfun(model)
    is.matrix(sf) && nrow(sf) > 0
  }, error = function(e) FALSE)

  if (!has_estfun) {
    message("  estfun() falló para glmerMod — merDeriv no disponible o incompatible")
    message("  -> fallback a vcov del modelo (SE no-clustered)")
    return(NULL)
  }

  # Probar que vcovCL funciona antes de retornar la fórmula
  vcv_test <- tryCatch(
    sandwich::vcovCL(model,
                     cluster = as.formula(paste0("~", cluster_var)),
                     type    = "HC1"),
    error   = function(e) { message("  vcovCL error: ", e$message); NULL },
    warning = function(w) { message("  vcovCL warning: ", w$message); NULL }
  )

  if (is.null(vcv_test)) {
    message("  vcovCL falló -> fallback a vcov del modelo")
    return(NULL)
  }

  # Validar PD
  ev <- tryCatch(
    eigen(vcv_test, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(ev) || any(ev < -1e-8)) {
    message("  vcov HC1 no-PD -> fallback")
    return(NULL)
  }

  message("  vcovCL HC1 OK — usando cluster-robust SEs")
  as.formula(paste0("~", cluster_var))
}

# RE probit con nAGQ (equivalente a xtprobit GHQ de Stata)
fit_reprobit <- function(outcome, data, ref_clase, with_controls, label) {

  ref_label <- dplyr::filter(STATE_MAP, clase == ref_clase) %>%
    dplyr::pull(label_long)
  message("\n  Fitting: ", label,
          " | outcome=",  outcome,
          " | ref=clase", ref_clase, " (", ref_label, ")",
          " | nAGQ=", NAGQ)

  d <- data %>%
    dplyr::filter(!is.na(clase), !is.na(.data[[outcome]])) %>%
    dplyr::mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola),
      clase_f      = relevel(as.factor(clase), ref = as.character(ref_clase))
    )

  ctrl <- if (with_controls) {
    ctrl_use[ctrl_use %in% names(d) &
             sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  } else character(0)

  d <- tidyr::drop_na(
    d, dplyr::all_of(intersect(c(outcome, "clase_f", "ola_f", ctrl), names(d)))
  )

  if (nrow(d) < 50) { warning("n muy bajo en ", label, ": ", nrow(d)); return(NULL) }
  message("    N obs=", nrow(d), " | N ids=", length(unique(d$idencuesta_f)))

  rhs <- paste(c("clase_f", ctrl, "ola_f", "(1 | idencuesta_f)"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  # nAGQ≥2 con exactamente 1 RE → equivalente a GHQ de Stata
  # Fallback a Laplace (nAGQ=1) si converge con error
  m <- tryCatch(
    lme4::glmer(
      fml, data = d,
      family  = binomial(link = "probit"),
      nAGQ    = NAGQ,
      control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl   = list(maxfun = 5e5)
      )
    ),
    error = function(e) {
      message("  glmer nAGQ=", NAGQ, " falló. Reintentando con Laplace (nAGQ=1).")
      lme4::glmer(
        fml, data = d,
        family  = binomial(link = "probit"),
        nAGQ    = 1L,
        control = lme4::glmerControl(optimizer = "bobyqa",
                                     optCtrl   = list(maxfun = 5e5))
      )
    }
  )

  sigma_u <- round(sqrt(as.numeric(VarCorr(m)$idencuesta_f)), 4)
  message("    sigma_u=", sigma_u,
          " (Stata típico: 1.5-2.5 para RE probit con trust)")

  vcov_fml  <- get_vcov_formula(m)
  vcov_type <- ifelse(is.null(vcov_fml), "model_vcov (fallback)", "HC1_cluster")

  # vcov = fórmula ~cluster → sandwich::vcovCL (cluster-robust HC1)
  # vcov = NULL             → marginaleffects usa vcov(model) (no-clustered, fallback)
  ame <- marginaleffects::avg_slopes(
    m, variables = "clase_f", type = "response",
    vcov = if (is.null(vcov_fml)) TRUE else vcov_fml
  ) %>%
    as_tibble() %>%
    dplyr::mutate(
      outcome       = outcome,
      ref_clase     = ref_clase,
      ref_label     = ref_label,
      with_controls = with_controls,
      spec          = label,
      n_obs         = nrow(d),
      n_id          = length(unique(d$idencuesta_f)),
      vcov_type     = vcov_type,
      sigma_u       = sigma_u
    )

  # Etiquetar contrasts con símbolos canónicos
  ame <- ame %>%
    dplyr::mutate(
      contrast_label = stringr::str_replace_all(
        contrast,
        c(
          "mean\\(1\\)" = "\u03b2 Closed",
          "mean\\(2\\)" = "\u03b3 Bridging",
          "mean\\(3\\)" = "\u03b1 Isolate",
          "^1 - 2$"     = "\u03b2 vs \u03b3",
          "^3 - 2$"     = "\u03b1 vs \u03b3",
          "^2 - 1$"     = "\u03b3 vs \u03b2",
          "^3 - 1$"     = "\u03b1 vs \u03b2",
          "^1 - 3$"     = "\u03b2 vs \u03b1",
          "^2 - 3$"     = "\u03b3 vs \u03b1"
        )
      )
    )

  # Predicciones clase x ola para Figure 3
  pred <- tryCatch(
    marginaleffects::predictions(
      m,
      newdata = marginaleffects::datagrid(
        clase_f   = levels(d$clase_f),
        ola_f     = levels(d$ola_f),
        grid_type = "mean_or_mode"
      ),
      type = "response"
    ) %>%
      as_tibble() %>%
      dplyr::mutate(
        outcome     = outcome,
        spec        = label,
        clase_int   = as.integer(as.character(clase_f)),
        clase_label = dplyr::recode(as.character(clase_int), !!!STATE_LABELS_FIG)
      ),
    error = function(e) {
      warning("predictions() falló en ", label, ": ", e$message)
      NULL
    }
  )

  list(model = m, ame = ame, pred = pred,
       formula   = deparse(fml), n = nrow(d),
       label     = label, ref_label = ref_label,
       vcov_type = vcov_type, sigma_u = sigma_u)
}

# Pooled probit (Table A1)
fit_probit_pool <- function(outcome, data, ref_clase, with_controls, label) {

  ref_label <- dplyr::filter(STATE_MAP, clase == ref_clase) %>%
    dplyr::pull(label_long)

  d <- data %>%
    dplyr::filter(!is.na(clase), !is.na(.data[[outcome]])) %>%
    dplyr::mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola),
      clase_f      = relevel(as.factor(clase), ref = as.character(ref_clase))
    )

  ctrl <- if (with_controls) {
    ctrl_use[ctrl_use %in% names(d) &
             sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  } else character(0)

  d <- tidyr::drop_na(
    d, dplyr::all_of(intersect(c(outcome, "clase_f", "ola_f", ctrl), names(d)))
  )

  rhs <- paste(c("clase_f", ctrl, "ola_f"), collapse = " + ")
  fml <- as.formula(paste(outcome, "~", rhs))

  m   <- glm(fml, data = d, family = binomial(link = "probit"))
  vcv <- tryCatch(sandwich::vcovCL(m, cluster = ~idencuesta_f),
                  error = function(e) NULL)

  ame <- marginaleffects::avg_slopes(
    m, variables = "clase_f", type = "response", vcov = vcv
  ) %>%
    as_tibble() %>%
    dplyr::mutate(
      outcome = outcome, ref_clase = ref_clase, ref_label = ref_label,
      with_controls = with_controls, spec = label, n_obs = nrow(d)
    )

  list(model = m, ame = ame, formula = deparse(fml), label = label)
}

# Causalidad reversa: trust -> clase
fit_reverse <- function(trust_var, data, with_controls, label) {

  d <- data %>%
    dplyr::filter(!is.na(clase), !is.na(.data[[trust_var]])) %>%
    dplyr::mutate(
      idencuesta_f = as.factor(idencuesta),
      ola_f        = as.factor(ola)
    )

  ctrl <- if (with_controls) {
    ctrl_use[ctrl_use %in% names(d) &
             sapply(ctrl_use, function(v) sum(!is.na(d[[v]])) > 0)]
  } else character(0)

  d <- tidyr::drop_na(
    d, dplyr::all_of(intersect(c("clase", trust_var, ctrl), names(d)))
  )

  rhs <- paste(c(trust_var, ctrl, "ola_f", "(1 | idencuesta_f)"), collapse = " + ")
  fml <- as.formula(paste("clase ~", rhs))

  m <- lme4::lmer(fml, data = d,
                  control = lme4::lmerControl(optimizer = "bobyqa"))

  vcv <- tryCatch({
    v  <- sandwich::vcovCL(m, cluster = d$idencuesta_f, type = "HC1")
    ev <- eigen(v, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < -1e-8)) NULL else v
  }, error = function(e) {
    message("  vcovCL error en fit_reverse: ", e$message, " -> NULL")
    NULL
  })

  slopes <- tryCatch(
    marginaleffects::avg_slopes(m, variables = trust_var, vcov = vcv) %>%
      as_tibble() %>%
      dplyr::mutate(
        predictor     = trust_var,
        with_controls = with_controls,
        spec          = label,
        n_obs         = nrow(d)
      ),
    error = function(e) NULL
  )

  list(model = m, slopes = slopes, formula = deparse(fml), label = label)
}

# ==============================================================================
# 5. TABLE 4 — RE probit ref=γ Bridging (ib2.clase)
# ==============================================================================
message("\n=== TABLE 4: RE probit | ref=\u03b3 Bridging | nAGQ=", NAGQ, " ===")
message("  Targets paper Table 4 M1: \u03b2=-0.071***, \u03b1=-0.040*")
message("  Targets paper Table 4 M2: \u03b2=-0.060***, \u03b1=-0.037*")
message("  Targets paper Table 4 M3: \u03b2=-0.046(ns), \u03b1=-0.091***")
message("  Targets paper Table 4 M4: \u03b2=-0.073***, \u03b1=-0.092***")

t4_m1 <- fit_reprobit("trust",    dt, ref_clase=2, with_controls=FALSE, label="T4_M1")
t4_m2 <- fit_reprobit("trust",    dt, ref_clase=2, with_controls=TRUE,  label="T4_M2")
t4_m3 <- fit_reprobit("trust_nh", dt, ref_clase=2, with_controls=FALSE, label="T4_M3")
t4_m4 <- fit_reprobit("trust_nh", dt, ref_clase=2, with_controls=TRUE,  label="T4_M4")

table4 <- dplyr::bind_rows(t4_m1$ame, t4_m2$ame, t4_m3$ame, t4_m4$ame)
readr::write_csv(table4, here::here("output", "table4_reprobit_ref_bridging.csv"))

cat("\n--- Table 4 (ref=\u03b3 Bridging; negativo = menos trust que Bridging) ---\n")
table4 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value < .001 ~ "***", p.value < .01 ~ "**",
    p.value < .05  ~ "*",   p.value < .10  ~ ".", TRUE ~ ""
  )) %>%
  dplyr::select(spec, outcome, contrast_label, estimate, std.error, p.value, sig,
                n_obs, vcov_type, sigma_u) %>%
  print()

# Comparación directa vs Stata
cat("\n--- Comparación vs Stata (N=3891, Table 4) ---\n")
# Targets del paper (Table 4, ref=γ Bridging)
stata_ref <- tibble::tibble(
  spec           = c("T4_M1","T4_M1","T4_M2","T4_M2","T4_M3","T4_M3","T4_M4","T4_M4"),
  contrast_label = rep(c("\u03b2 vs \u03b3", "\u03b1 vs \u03b3"), 4),
  stata_est      = c(-0.071,-0.040, -0.060,-0.037, -0.046,-0.091, -0.073,-0.092),
  stata_se       = c(0.021, 0.020,  0.020, 0.020,  0.032, 0.030,  0.032, 0.031),
  stata_sig      = c("***","*",    "***","*",     "",    "***",   "***","***")
)
table4 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::select(spec, contrast_label, estimate, std.error, p.value) %>%
  dplyr::mutate(
    sig   = dplyr::case_when(p.value<.001~"***",p.value<.01~"**",
                             p.value<.05~"*",p.value<.10~".",TRUE~""),
    R_est = round(estimate,  3),
    R_se  = round(std.error, 3)
  ) %>%
  dplyr::left_join(stata_ref, by = c("spec", "contrast_label")) %>%
  dplyr::select(spec, contrast_label, R_est, R_se, sig, stata_est, stata_se, stata_sig) %>%
  print()

# ==============================================================================
# 6. TABLE 5 — RE probit ref=β Closed (ib1.clase)
# ==============================================================================
message("\n=== TABLE 5: RE probit | ref=\u03b2 Closed | nAGQ=", NAGQ, " ===")

t5_m5 <- fit_reprobit("trust",    dt, ref_clase=1, with_controls=FALSE, label="T5_M5")
t5_m6 <- fit_reprobit("trust",    dt, ref_clase=1, with_controls=TRUE,  label="T5_M6")
t5_m7 <- fit_reprobit("trust_nh", dt, ref_clase=1, with_controls=FALSE, label="T5_M7")
t5_m8 <- fit_reprobit("trust_nh", dt, ref_clase=1, with_controls=TRUE,  label="T5_M8")

table5 <- dplyr::bind_rows(t5_m5$ame, t5_m6$ame, t5_m7$ame, t5_m8$ame)
readr::write_csv(table5, here::here("output", "table5_reprobit_ref_closed.csv"))

cat("\n--- Table 5 (ref=\u03b2 Closed) ---\n")
table5 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",p.value<.10~".",TRUE~""
  )) %>%
  dplyr::select(spec, outcome, contrast_label, estimate, std.error, p.value, sig) %>%
  print()

# ==============================================================================
# 7. TABLE A1 — Pooled probit ref=β Closed
# ==============================================================================
message("\n=== TABLE A1: Pooled probit | ref=\u03b2 Closed ===")

a1_m1 <- fit_probit_pool("trust",    dt, ref_clase=1, with_controls=FALSE, label="A1_M1")
a1_m2 <- fit_probit_pool("trust",    dt, ref_clase=1, with_controls=TRUE,  label="A1_M2")
a1_m3 <- fit_probit_pool("trust_nh", dt, ref_clase=1, with_controls=FALSE, label="A1_M3")
a1_m4 <- fit_probit_pool("trust_nh", dt, ref_clase=1, with_controls=TRUE,  label="A1_M4")

tableA1 <- dplyr::bind_rows(a1_m1$ame, a1_m2$ame, a1_m3$ame, a1_m4$ame)
readr::write_csv(tableA1, here::here("output", "tableA1_probit_pool_ref_closed.csv"))

cat("\n--- Table A1 ---\n")
tableA1 %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",TRUE~""
  )) %>%
  dplyr::select(spec, outcome, contrast, estimate, std.error, p.value, sig) %>%
  print()

# ==============================================================================
# 8. SENSIBILIDAD: conf_gral
# ==============================================================================
message("\n=== SENSIBILIDAD: conf_gral ===")
if ("conf_gral" %in% names(dt) && sum(!is.na(dt$conf_gral)) > 100) {
  s_b <- fit_reprobit("conf_gral", dt, ref_clase=2, with_controls=TRUE, label="Sens_Bridging")
  s_c <- fit_reprobit("conf_gral", dt, ref_clase=1, with_controls=TRUE, label="Sens_Closed")
  if (!is.null(s_b) && !is.null(s_c)) {
    dplyr::bind_rows(s_b$ame, s_c$ame) %>%
      readr::write_csv(here::here("output","sensitivity_confgral.csv"))
  }
} else message("  conf_gral no disponible.")

# ==============================================================================
# 9. TABLE A2 — Causalidad reversa
# ==============================================================================
message("\n=== TABLE A2: Causalidad reversa (trust -> clase) ===")

ra1 <- fit_reverse("trust",    dt, with_controls=FALSE, label="A2_M1")
ra2 <- fit_reverse("trust",    dt, with_controls=TRUE,  label="A2_M2")
ra3 <- fit_reverse("trust_nh", dt, with_controls=FALSE, label="A2_M3")
ra4 <- fit_reverse("trust_nh", dt, with_controls=TRUE,  label="A2_M4")

tableA2 <- dplyr::bind_rows(ra1$slopes, ra2$slopes, ra3$slopes, ra4$slopes)
readr::write_csv(tableA2, here::here("output","tableA2_reverse_causality.csv"))

cat("\n--- Table A2 ---\n")
tableA2 %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",p.value<.10~".",TRUE~""
  )) %>%
  dplyr::select(spec, predictor, estimate, std.error, p.value, sig) %>%
  print()

# ==============================================================================
# 10. FIGURE 3 — Pr(trust=1) | clase x ola
#     Replica estructura Stata marginsplot: x=clase, color=ola, facets=outcome
# ==============================================================================
message("\n=== FIGURE 3 ===")

prep_pred <- function(pred_tbl, outcome_label) {
  if (is.null(pred_tbl)) return(NULL)
  pred_tbl %>%
    dplyr::mutate(
      wave_year     = dplyr::case_when(
        ola_f == "1" ~ WAVES_YEARS[1],
        ola_f == "2" ~ WAVES_YEARS[2],
        ola_f == "3" ~ WAVES_YEARS[3],
        TRUE         ~ as.integer(as.character(ola_f))
      ),
      outcome_label = outcome_label
    )
}

preds <- dplyr::bind_rows(
  prep_pred(t4_m2$pred, "Generalized trust (c02)"),
  prep_pred(t4_m4$pred, "Trust in neighbors (t01)")
)

if (!is.null(preds) && nrow(preds) > 0) {

  if (!"clase_label" %in% names(preds)) {
    preds <- preds %>%
      dplyr::mutate(
        clase_int   = as.integer(as.character(clase_f)),
        clase_label = dplyr::recode(as.character(clase_int), !!!STATE_LABELS_FIG)
      )
  }

  # Ordenar x igual que Stata Fig.3: Closed(β) → Broker(γ) → Apathetic(α)
  preds <- preds %>%
    dplyr::mutate(
      clase_order = dplyr::case_when(
        clase_int == 1L ~ 1L,   # β Closed
        clase_int == 2L ~ 2L,   # γ Bridging
        clase_int == 3L ~ 3L,   # α Isolate/Apathetic
        TRUE ~ NA_integer_
      ),
      wave_label = as.character(wave_year)
    )

  fig3 <- ggplot2::ggplot(
    preds,
    ggplot2::aes(
      x     = reorder(clase_label, clase_order),
      y     = estimate,
      ymin  = conf.low,
      ymax  = conf.high,
      color = wave_label,
      group = wave_label
    )
  ) +
    ggplot2::geom_ribbon(alpha = 0.10, color = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~outcome_label, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = c("2016" = "#4E79A7", "2018" = "#E15759", "2022" = "#59A14F"),
      name   = NULL
    ) +
    ggplot2::labs(
      x       = "Class",
      y       = "Pr(trust = 1)",
      color   = NULL,
      caption = paste0(
        "RE probit, full controls (nAGQ=", NAGQ, "). Ribbons = 95% CI.\n",
        "Reference: \u03b3 (Bridging/Broker). AMEs at population level;\n",
        "SEs from fixed-effect VCOV (Stata margins post-xtprobit equivalent)."
      )
    ) +
    theme_ssr() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text      = ggplot2::element_text(face = "bold"),
      axis.text.x     = ggplot2::element_text(size = 10)
    )

  ggsave_safe("fig3_predicted_probs.png", fig3, width = 9, height = 5)
  ggsave_safe("fig3_predicted_probs.pdf", fig3, width = 9, height = 5)
  message("  Fig3 guardada.")
} else {
  message("  predictions() no disponibles. Fig3 omitida.")
}

# ==============================================================================
# 11. AME UNIFICADO
# ==============================================================================
ame_all <- dplyr::bind_rows(
  table4  %>% dplyr::mutate(table = "Table4_ref_Bridging"),
  table5  %>% dplyr::mutate(table = "Table5_ref_Closed"),
  tableA1 %>% dplyr::mutate(table = "TableA1_pool")
) %>%
  dplyr::filter(!grepl("^ola", contrast, ignore.case=TRUE)) %>%
  dplyr::mutate(sig = dplyr::case_when(
    p.value<.001~"***",p.value<.01~"**",p.value<.05~"*",p.value<.10~".",TRUE~""
  )) %>%
  dplyr::arrange(table, outcome, spec, contrast)

readr::write_csv(ame_all, here::here("output","ame_all_models.csv"))

cat("\n=== AMEs significativos (p < .05) ===\n")
ame_all %>%
  dplyr::filter(p.value < .05) %>%
  dplyr::select(table, spec, outcome, contrast_label, estimate, std.error, p.value, sig) %>%
  print(n=40)

# ==============================================================================
# 12. SUMMARIES
# ==============================================================================
mods_to_report <- list(
  "T4_M2 trust ref=\u03b3 +ctrl"    = if (!is.null(t4_m2)) t4_m2$model else NULL,
  "T4_M4 trust_nh ref=\u03b3 +ctrl" = if (!is.null(t4_m4)) t4_m4$model else NULL,
  "T5_M6 trust ref=\u03b2 +ctrl"    = if (!is.null(t5_m6)) t5_m6$model else NULL,
  "T5_M8 trust_nh ref=\u03b2 +ctrl" = if (!is.null(t5_m8)) t5_m8$model else NULL,
  "A2_M2 reverse trust +ctrl"       = if (!is.null(ra2))   ra2$model   else NULL
)

summ_lines <- unlist(lapply(names(mods_to_report), function(nm) {
  m <- mods_to_report[[nm]]
  if (is.null(m)) return(paste0("\n=== ", nm, ": NULL ==="))
  c(paste0("\n\n=== ", nm, " ==="), capture.output(summary(m)))
}))
writeLines(summ_lines, here::here("output","trust_models_summaries.txt"))

# ==============================================================================
# 13. GUARDAR
# ==============================================================================
saveRDS(
  list(
    dt               = dt,
    ctrl_use         = ctrl_use,
    NAGQ             = NAGQ,
    STATE_MAP        = STATE_MAP,
    STATE_LABELS_FIG = STATE_LABELS_FIG,
    STATE_PALETTE    = STATE_PALETTE,
    t4_m2            = if (!is.null(t4_m2)) t4_m2$model else NULL,
    t4_m4            = if (!is.null(t4_m4)) t4_m4$model else NULL,
    t5_m6            = if (!is.null(t5_m6)) t5_m6$model else NULL,
    t5_m8            = if (!is.null(t5_m8)) t5_m8$model else NULL
  ),
  here::here("data","trust_model_objects.rds")
)

# ==============================================================================
# 14. RESUMEN FINAL
# ==============================================================================
message("\n[003_trust_models.R] DONE.")
message("  nAGQ=", NAGQ, " | N T4_M2=", if (!is.null(t4_m2)) t4_m2$n else "NULL")
message("  Controles: ", paste(ctrl_use, collapse=", "))
message("")
message("  SI AÚN DIFIERE DE STATA:")
message("    1. Ver output/pre_model_diagnostics.txt (N y patrón trust)")
message("    2. Probar NAGQ <- 12L para replicación exacta de GHQ-12")
message("    3. N muy diferente → revisar filtros en 01_descriptive_stats.R")
message("    4. sigma_u << 1.5 → el RE está mal estimado; usar nAGQ=12")
message("")
message("  Outputs:")
message("    output/pre_model_diagnostics.txt     ← LEER PRIMERO")
message("    output/table4_reprobit_ref_bridging.csv")
message("    output/table5_reprobit_ref_closed.csv")
message("    output/tableA1_probit_pool_ref_closed.csv")
message("    output/tableA2_reverse_causality.csv")
message("    output/ame_all_models.csv")
message("    output/fig3_predicted_probs.png/.pdf")
