# ==============================================================================
# check_coding.R — Verifica distribución RAW de c12 en ELSOC original
# Corre desde la raíz del proyecto: Rscript check_coding.R
# ==============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(here) })

ITEMS   <- c(paste0("c12_0", 1:8), "c12_09")
DOMAINS <- c("nhg","religious","political","union",
              "professional","charity","sport","student","other")
NA_CODES <- c(-999L, -888L, -777L, -666L)

# ------------------------------------------------------------------------------
# 1. Cargar ELSOC
# ------------------------------------------------------------------------------
env  <- new.env()
nms  <- load(here("data/ELSOC_Long.RData"), envir=env)
elsoc <- get(if ("elsoc_long_2016_2022" %in% nms) "elsoc_long_2016_2022" else nms[1],
             envir=env)

dt <- elsoc %>%
  filter(ola %in% c(1,3,6)) %>%
  { if ("muestra"       %in% names(.)) filter(., muestra==1)       else . } %>%
  { if ("tipo_atricion" %in% names(.)) filter(., tipo_atricion==1) else . } %>%
  mutate(across(all_of(ITEMS), ~ replace(.x, .x %in% NA_CODES, NA)))

cat("Obs totales (olas 1/3/6, filtros aplicados):", nrow(dt),
    "| N personas:", n_distinct(dt$idencuesta), "\n\n")

# ------------------------------------------------------------------------------
# 2. Distribución RAW de c12 por dominio
# ------------------------------------------------------------------------------
cat("=== DISTRIBUCIÓN RAW c12 POR DOMINIO ===\n")
cat("  1 = No es miembro\n  2 = Miembro INACTIVO\n  3 = Miembro ACTIVO\n\n")

raw_tbl <- lapply(seq_along(ITEMS), function(i) {
  x <- dt[[ITEMS[i]]]
  n_total <- sum(!is.na(x))
  tibble(
    dominio    = DOMAINS[i],
    `1_no`     = sum(x==1L, na.rm=TRUE),
    `2_inactivo`= sum(x==2L, na.rm=TRUE),
    `3_activo`  = sum(x==3L, na.rm=TRUE),
    NA_         = sum(is.na(x)),
    pct_inactivo= round(100*sum(x==2L,na.rm=TRUE)/n_total, 1),
    pct_activo  = round(100*sum(x==3L,na.rm=TRUE)/n_total, 1),
    pct_member_any  = round(100*sum(x>=2L,na.rm=TRUE)/n_total, 1),  # c12>=2
    pct_member_ao   = round(100*sum(x==3L,na.rm=TRUE)/n_total, 1)   # c12==3
  )
}) %>% bind_rows()

print(as.data.frame(raw_tbl))

# ------------------------------------------------------------------------------
# 3. Qué CAMBIA entre any_member y active_only: los c12==2
# ------------------------------------------------------------------------------
cat("\n=== IMPACTO DEL CAMBIO DE CODIFICACIÓN ===\n")
cat("Solo cambian las personas con c12 == 2 (inactivo):\n")
cat("  any_member: c12==2 → 1 (cuenta como miembro)\n")
cat("  active_only: c12==2 → 0 (NO cuenta como miembro)\n\n")

impact_tbl <- lapply(seq_along(ITEMS), function(i) {
  x <- dt[[ITEMS[i]]]
  tibble(
    dominio = DOMAINS[i],
    n_c12_2  = sum(x==2L, na.rm=TRUE),
    pct_c12_2= round(100*sum(x==2L,na.rm=TRUE)/sum(!is.na(x)),1)
  )
}) %>% bind_rows()

print(as.data.frame(impact_tbl))
cat("\nTotal c12==2 (obs que cambian de 1→0):",
    sum(sapply(ITEMS, function(i) sum(dt[[i]]==2L, na.rm=TRUE))),
    "de",
    sum(sapply(ITEMS, function(i) sum(!is.na(dt[[i]]))), na.rm=TRUE),
    "obs totales\n")

# ------------------------------------------------------------------------------
# 4. A nivel de PERSONA: ¿cuántos dominios activos vs. any_member?
# ------------------------------------------------------------------------------
cat("\n=== A NIVEL DE PERSONA-OLA: n_domains bajo cada codificación ===\n")
dt_per <- dt %>%
  mutate(
    n_any = rowSums(across(all_of(ITEMS), ~as.integer(!is.na(.x) & .x>=2L)), na.rm=TRUE),
    n_ao  = rowSums(across(all_of(ITEMS), ~as.integer(!is.na(.x) & .x==3L)), na.rm=TRUE),
    diff  = n_any - n_ao   # cuántos dominios se "pierden" al pasar a active_only
  )

cat("Distribución n_domains ANY_MEMBER (c12>=2):\n")
print(table(dt_per$n_any))
cat("\nDistribución n_domains ACTIVE_ONLY (c12==3):\n")
print(table(dt_per$n_ao))
cat("\nDiferencia (n_any - n_ao) — dominios que se pierden por cambio de codificación:\n")
print(table(dt_per$diff))

cat("\n% de personas-ola con al menos 1 dominio que cambia (c12==2 en algún item):\n")
cat(round(100*mean(dt_per$diff > 0), 1), "%\n")

cat("\n=== VERIFICACIÓN CRUZADA ===\n")
cat("Pregunta: ¿puede alguien tener n_ao > n_any? (NO debería, c12==3 es subconjunto de c12>=2)\n")
cat("Casos con n_ao > n_any:", sum(dt_per$n_ao > dt_per$n_any), "\n")

cat("\n[check_coding.R DONE]\n")
