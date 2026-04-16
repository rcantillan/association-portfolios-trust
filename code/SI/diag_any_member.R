# ==============================================================================
# diag_any_member.R — Diagnóstico pipeline con codificación ANY_MEMBER (c12 >= 2)
# ------------------------------------------------------------------------------
# Corre de forma independiente para verificar manualmente:
#   (1) Distribución raw de c12 por dominio
#   (2) Recodificación binaria c12 >= 2
#   (3) Conteo de dominios por persona
#   (4) Clasificación LMM (dt_states_cov.rds)
#   (5) Trust por clase
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(here); library(readr)
})

CODING <- "any_member"
cat("=================================================================\n")
cat("DIAGNÓSTICO CODIFICACIÓN:", CODING, "\n")
cat("  Regla: c12 >= 2 (miembro inactivo o activo = 1; no miembro = 0)\n")
cat("=================================================================\n\n")

items    <- paste0("c12_0", 1:8)
dominios <- c("nhg","religious","political","union",
              "professional","charity","sport","student")

# ------------------------------------------------------------------------------
# 1. DISTRIBUCIÓN RAW DE c12 POR DOMINIO (ELSOC Long)
# ------------------------------------------------------------------------------
stopifnot(file.exists(here("data/ELSOC_Long.RData")))
obj_names <- load(here("data/ELSOC_Long.RData"))
a <- get(if ("elsoc_long_2016_2022" %in% obj_names) "elsoc_long_2016_2022" else obj_names[1])
a <- a %>% filter(ola %in% c(1,3,6))
if ("muestra"      %in% names(a)) a <- a %>% filter(muestra == 1)
if ("tipo_atricion" %in% names(a)) a <- a %>% filter(tipo_atricion == 1)

# Reemplazar códigos de missing
a <- a %>% mutate(across(all_of(items),
                         ~ replace(.x, .x %in% c(-999,-888,-666), NA)))

cat("--- 1. DISTRIBUCIÓN RAW c12 (1=no miembro, 2=inactivo, 3=activo) ---\n")
raw_dist <- a %>%
  select(all_of(items)) %>%
  pivot_longer(everything(), names_to="item", values_to="valor") %>%
  count(item, valor) %>%
  mutate(
    dominio = dominios[match(item, items)],
    pct = round(n / sum(n[!is.na(valor)]) * 100, 1)
  ) %>%
  filter(!is.na(valor)) %>%
  pivot_wider(id_cols=c(item,dominio), names_from=valor,
              values_from=n, names_prefix="c12_") %>%
  mutate(total = rowSums(across(starts_with("c12_")), na.rm=TRUE),
         pct_member  = round(100 * (c12_2 + c12_3) / total, 1),
         pct_active  = round(100 * c12_3 / total, 1),
         pct_inactive= round(100 * c12_2 / total, 1))

print(raw_dist %>% select(dominio, c12_1, c12_2, c12_3,
                           pct_member, pct_inactive, pct_active))

# ------------------------------------------------------------------------------
# 2. RECODIFICACIÓN BINARIA: c12 >= 2 → 1
# ------------------------------------------------------------------------------
cat("\n--- 2. RECODIFICACIÓN any_member: c12 >= 2 → 1 ---\n")
a_bin <- a %>%
  mutate(across(all_of(items), ~ as.integer(!is.na(.x) & .x >= 2)))

cat("Distribución de SUMAS por persona (n_dominios con membresía):\n")
a_bin <- a_bin %>%
  mutate(n_member_any = rowSums(across(all_of(items)), na.rm=TRUE))
print(table(a_bin$n_member_any, useNA="ifany"))

cat("\nProp con 0, 1-2, 3-4, 5+ dominios:\n")
a_bin %>%
  mutate(grupo = case_when(
    n_member_any == 0   ~ "0 dominios",
    n_member_any <= 2   ~ "1-2 dominios",
    n_member_any <= 4   ~ "3-4 dominios",
    TRUE                ~ "5+ dominios"
  )) %>%
  count(grupo) %>%
  mutate(pct = round(100*n/sum(n),1)) %>%
  print()

# ------------------------------------------------------------------------------
# 3. PERSONAS EN PANEL BALANCEADO (mismo N que el LMM)
# ------------------------------------------------------------------------------
cat("\n--- 3. PANEL BALANCEADO (3 olas) ---\n")
ids_3olas <- a_bin %>%
  filter(ola %in% c(1,3,6)) %>%
  count(idencuesta) %>%
  filter(n == 3) %>%
  pull(idencuesta)

dt_bal <- a_bin %>%
  filter(idencuesta %in% ids_3olas) %>%
  select(idencuesta, ola, all_of(items), n_member_any)

cat("N personas:", length(ids_3olas), "| N obs:", nrow(dt_bal), "\n")

cat("\nDistribución n_dominios (any_member) en panel balanceado:\n")
print(table(dt_bal$n_member_any))

# ------------------------------------------------------------------------------
# 4. CLASIFICACIÓN LMM (cargada desde disco)
# ------------------------------------------------------------------------------
cat("\n--- 4. CLASIFICACIÓN LMM (any_member: dt_states_cov.rds) ---\n")
dt_lmm <- readRDS(here("data/dt_states_cov.rds")) %>%
  mutate(idencuesta = as.character(idencuesta))

cat("Tamaño de clases:\n")
print(dt_lmm %>% count(position) %>% mutate(pct=round(100*n/sum(n),1)))

cat("\nn_member_any MEDIO por clase (membresías any_member):\n")
dt_lmm_bin <- dt_lmm %>%
  mutate(idencuesta = as.character(idencuesta)) %>%
  left_join(
    dt_bal %>% mutate(idencuesta=as.character(idencuesta)),
    by = c("idencuesta","ola")
  )
print(dt_lmm_bin %>%
  group_by(position) %>%
  summarise(
    mean_n_any  = round(mean(n_member_any, na.rm=TRUE), 2),
    median_n    = median(n_member_any, na.rm=TRUE),
    pct_0 = round(mean(n_member_any==0, na.rm=TRUE)*100, 1),
    pct_1 = round(mean(n_member_any==1, na.rm=TRUE)*100, 1),
    pct_2plus = round(mean(n_member_any>=2, na.rm=TRUE)*100, 1)
  ))

cat("\nProporción observada c12>=2 por clase y dominio:\n")
print(dt_lmm_bin %>%
  group_by(position) %>%
  summarise(across(all_of(items), ~ round(mean(.x, na.rm=TRUE), 3))))

# ------------------------------------------------------------------------------
# 5. TRUST POR CLASE
# ------------------------------------------------------------------------------
cat("\n--- 5. TRUST POR CLASE (any_member) ---\n")
dt_anal <- readRDS(here("data/dt_analysis.rds")) %>%
  rename(idencuesta=id) %>%
  mutate(idencuesta=as.character(idencuesta)) %>%
  select(idencuesta, ola, trust, trust_nh)

dt_trust <- dt_lmm %>%
  left_join(dt_anal, by=c("idencuesta","ola"))

print(dt_trust %>%
  group_by(position) %>%
  summarise(
    n = n(),
    trust_gen = round(mean(trust,    na.rm=TRUE)*100, 1),
    trust_nh  = round(mean(trust_nh, na.rm=TRUE)*100, 1)
  ))

cat("\n=== FIN diag_any_member.R ===\n")
