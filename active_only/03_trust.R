# ==============================================================================
# active_only / 03_trust.R
# Modelos RE probit de confianza por clase LMM (codificación ACTIVE_ONLY)
# Input:  output/dt_panel.rds  (trust, controles)
#         output/dt_states.rds (clases α/β/γ)
# Output: output/trust_ames.csv
#         output/trust_summary.txt
# ==============================================================================

options(marginaleffects_safe = FALSE)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr); library(here)
  library(lme4); library(marginaleffects); library(clubSandwich); library(sandwich)
})
if (!requireNamespace("merDeriv", quietly = TRUE))
  install.packages("merDeriv", quiet = TRUE)
suppressPackageStartupMessages(library(merDeriv))

ROOT    <- here::here()
DIR_OUT <- file.path(ROOT, "active_only", "output")

NAGQ     <- 12L
CONTROLS <- c("age","woman","education","swb","employed","couple")

cat("=== 03_trust.R | active_only | nAGQ =", NAGQ, "===\n\n")

# ------------------------------------------------------------------------------
# 1. Cargar y unir datos
# ------------------------------------------------------------------------------
dt_panel  <- readRDS(file.path(DIR_OUT, "dt_panel.rds")) %>%
  mutate(idencuesta = as.character(idencuesta))
dt_states <- readRDS(file.path(DIR_OUT, "dt_states.rds")) %>%
  mutate(idencuesta = as.character(idencuesta)) %>%
  select(idencuesta, ola, position, p_alpha, p_beta, p_gamma, p_max)

# Usar solo columnas de panel que no estén ya en states
cols_extra <- setdiff(names(dt_panel),
                      c("idencuesta","ola",
                        intersect(names(dt_panel), names(dt_states))))

dt <- left_join(dt_states,
                dt_panel %>% select(idencuesta, ola, all_of(cols_extra)),
                by = c("idencuesta","ola")) %>%
  mutate(
    # γ = referencia (clase 2 en factor)
    clase = case_when(position=="alpha"~3L, position=="beta"~1L, position=="gamma"~2L),
    clase_f = factor(clase, levels=c(2L,1L,3L)),   # ref = γ
    ola_f   = factor(ola)
  )

cat("N obs:", nrow(dt), "| N personas:", n_distinct(dt$idencuesta), "\n")
cat("Clases:\n")
print(dt %>% count(position) %>% mutate(pct=round(100*n/sum(n),1)))

cat("\nTrust por clase (descriptivo):\n")
dt_trust_desc <- left_join(
  dt_states %>% select(idencuesta, ola, position),
  dt_panel  %>% select(idencuesta, ola, trust, trust_nh),
  by=c("idencuesta","ola")
)
print(dt_trust_desc %>%
  group_by(position) %>%
  summarise(trust_gen=round(mean(trust,na.rm=TRUE)*100,1),
            trust_nh =round(mean(trust_nh,na.rm=TRUE)*100,1),
            n=n()))

# Estandarizar variables continuas (evita eigenvalue grande en glmer)
dt <- dt %>%
  mutate(age_z   = as.numeric(scale(age)),
         tinst_z = as.numeric(scale(tinst)))

# Reemplazar en CONTROLS: age → age_z, tinst → tinst_z
ctrl_ok <- CONTROLS[CONTROLS %in% names(dt) &
                    sapply(CONTROLS, function(v) sum(!is.na(dt[[v]])) > 0)]
ctrl_ok <- gsub("^age$",   "age_z",   ctrl_ok)
ctrl_ok <- gsub("^tinst$", "tinst_z", ctrl_ok)
ctrl_ok <- ctrl_ok[ctrl_ok %in% names(dt)]
cat("\nControles disponibles (age y tinst estandarizados):", paste(ctrl_ok, collapse=", "), "\n\n")

# ------------------------------------------------------------------------------
# 2. Función para estimar y extraer AMEs
# ------------------------------------------------------------------------------
fit_model <- function(outcome, controls, label) {
  d <- dt %>%
    drop_na(all_of(c(outcome, "clase_f", controls))) %>%
    mutate(idencuesta_f = as.factor(idencuesta))

  fml <- as.formula(
    if (length(controls) == 0)
      paste0(outcome, " ~ clase_f + ola_f + (1|idencuesta_f)")
    else
      paste0(outcome, " ~ clase_f + ", paste(controls, collapse="+"),
             " + ola_f + (1|idencuesta_f)")
  )

  cat("  Fitting", label, "| N =", nrow(d), "\n")

  tryCatch({
    m <- tryCatch(
      glmer(fml, data=d, family=binomial("probit"), nAGQ=NAGQ,
            control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5))),
      error = function(e) {
        cat("    ⚠ nAGQ=", NAGQ, "falló; reintentando con Laplace (nAGQ=1)\n")
        glmer(fml, data=d, family=binomial("probit"), nAGQ=1L,
              control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5)))
      }
    )

    # Convergencia
    conv <- m@optinfo$conv$lme4
    if (length(conv) > 0) cat("    ⚠ Convergencia:", paste(unlist(conv$messages), collapse=" | "), "\n")

    # sigma_u (debe ser ~1.5-2.5 para RE probit bien estimado)
    sigma_u <- round(sqrt(as.numeric(VarCorr(m)$idencuesta_f)), 4)
    cat("    sigma_u =", sigma_u, "\n")

    vcv <- tryCatch({
      v <- as.matrix(clubSandwich::vcovCR(m, cluster=d$idencuesta_f, type="CR1"))
      ev <- eigen(v, symmetric=TRUE, only.values=TRUE)$values
      if (any(ev < -1e-8)) NULL else v
    }, error   = function(e) { cat("    ⚠ vcovCR error:", e$message, "\n"); NULL },
       warning = function(w) { cat("    ⚠ vcovCR warning:", w$message, "\n"); NULL })

    if (is.null(vcv)) cat("    ⚠ Cluster-robust vcov falló; usando model-based\n")

    ame <- avg_slopes(m, variables="clase_f", type="response",
                      vcov=if (!is.null(vcv)) vcv else TRUE)

    ame %>% as_tibble() %>%
      mutate(modelo=label, outcome=outcome,
             controls=if(length(controls)>0) "yes" else "no",
             N=nrow(d)) %>%
      select(modelo, outcome, controls, N, contrast,
             estimate, std.error, statistic, p.value, conf.low, conf.high)
  }, error=function(e) {
    warning("Error en ", label, ": ", e$message)
    tibble()
  })
}

# ------------------------------------------------------------------------------
# 3. Correr los 4 modelos
# ------------------------------------------------------------------------------
cat("--- Modelos de confianza generalizada ---\n")
m1 <- fit_model("trust",    character(0), "M1_gen_noctrl")
m2 <- fit_model("trust",    ctrl_ok,      "M2_gen_ctrl")

cat("\n--- Modelos de confianza vecinal ---\n")
m3 <- fit_model("trust_nh", character(0), "M3_nh_noctrl")
m4 <- fit_model("trust_nh", ctrl_ok,      "M4_nh_ctrl")

results <- bind_rows(m1, m2, m3, m4) %>%
  mutate(sig = case_when(
    p.value < .001 ~ "***", p.value < .01 ~ "**",
    p.value < .05  ~ "*",   p.value < .10 ~ ".",
    TRUE           ~ ""
  ))

# ------------------------------------------------------------------------------
# 4. Mostrar y guardar resultados
# ------------------------------------------------------------------------------
cat("\n=== AMEs: clase vs γ (referencia) ===\n")
cat("  Contraste '1 - 2' = β vs γ\n")
cat("  Contraste '3 - 2' = α vs γ\n\n")
print(results %>% select(modelo, outcome, contrast, estimate, std.error, p.value, sig))

write_csv(results, file.path(DIR_OUT, "trust_ames.csv"))
cat("\n✓ Guardado:", file.path(DIR_OUT, "trust_ames.csv"), "\n")

# Resumen compacto
cat("\n=== RESUMEN COMPACTO (con controles) ===\n")
print(results %>%
  filter(controls=="yes") %>%
  select(outcome, contrast, estimate, std.error, p.value, sig) %>%
  mutate(across(c(estimate, std.error), ~round(.x,3)),
         p.value=round(p.value,3)))

cat("\n[03_trust.R DONE]\n")
