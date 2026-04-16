# ==============================================================================
# diag_active_only.R — Diagnóstico pipeline con codificación ACTIVE_ONLY (c12 == 3)
# ------------------------------------------------------------------------------
# Corre de forma independiente para verificar manualmente:
#   (1) Distribución raw de c12 por dominio
#   (2) Recodificación binaria c12 == 3
#   (3) Conteo de dominios por persona
#   (4) Clasificación LMM (dt_states_cov_active_only.rds)
#   (5) Trust por clase
#   (6) COMPARACIÓN DIRECTA con any_member: ¿quiénes cambian de clase?
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(here); library(readr)
})

CODING <- "active_only"
cat("=================================================================\n")
cat("DIAGNÓSTICO CODIFICACIÓN:", CODING, "\n")
cat("  Regla: c12 == 3 (miembro ACTIVO = 1; inactivo o no miembro = 0)\n")
cat("=================================================================\n\n")

items    <- paste0("c12_0", 1:8)
dominios <- c("nhg","religious","political","union",
              "professional","charity","sport","student")

# ------------------------------------------------------------------------------
# 1. DISTRIBUCIÓN RAW DE c12 POR DOMINIO
# ------------------------------------------------------------------------------
stopifnot(file.exists(here("data/ELSOC_Long.RData")))
obj_names <- load(here("data/ELSOC_Long.RData"))
a <- get(if ("elsoc_long_2016_2022" %in% obj_names) "elsoc_long_2016_2022" else obj_names[1])
a <- a %>% filter(ola %in% c(1,3,6))
if ("muestra"      %in% names(a)) a <- a %>% filter(muestra == 1)
if ("tipo_atricion" %in% names(a)) a <- a %>% filter(tipo_atricion == 1)
a <- a %>% mutate(across(all_of(items),
                         ~ replace(.x, .x %in% c(-999,-888,-666), NA)))

cat("--- 1. DISTRIBUCIÓN RAW c12 (idéntica a any_member — mismos datos) ---\n")
raw_dist <- a %>%
  select(all_of(items)) %>%
  pivot_longer(everything(), names_to="item", values_to="valor") %>%
  count(item, valor) %>%
  filter(!is.na(valor)) %>%
  mutate(dominio = dominios[match(item, items)]) %>%
  pivot_wider(id_cols=c(item,dominio), names_from=valor,
              values_from=n, names_prefix="c12_") %>%
  mutate(total = rowSums(across(starts_with("c12_")), na.rm=TRUE),
         pct_active  = round(100 * c12_3 / total, 1),
         pct_inactive= round(100 * c12_2 / total, 1),
         pct_none    = round(100 * c12_1 / total, 1))

print(raw_dist %>% select(dominio, c12_1, c12_2, c12_3, pct_none, pct_inactive, pct_active))

# ------------------------------------------------------------------------------
# 2. RECODIFICACIÓN BINARIA: c12 == 3 → 1
# ------------------------------------------------------------------------------
cat("\n--- 2. RECODIFICACIÓN active_only: c12 == 3 → 1 ---\n")
a_bin_ao <- a %>%
  mutate(across(all_of(items), ~ as.integer(!is.na(.x) & .x == 3)))

cat("Distribución de SUMAS por persona (n_dominios con membresía ACTIVA):\n")
a_bin_ao <- a_bin_ao %>%
  mutate(n_member_ao = rowSums(across(all_of(items)), na.rm=TRUE))
print(table(a_bin_ao$n_member_ao, useNA="ifany"))

cat("\nProp con 0, 1-2, 3-4, 5+ dominios activos:\n")
a_bin_ao %>%
  mutate(grupo = case_when(
    n_member_ao == 0   ~ "0 dominios",
    n_member_ao <= 2   ~ "1-2 dominios",
    n_member_ao <= 4   ~ "3-4 dominios",
    TRUE               ~ "5+ dominios"
  )) %>%
  count(grupo) %>%
  mutate(pct = round(100*n/sum(n),1)) %>%
  print()

# ------------------------------------------------------------------------------
# 3. PANEL BALANCEADO
# ------------------------------------------------------------------------------
cat("\n--- 3. PANEL BALANCEADO (3 olas) ---\n")
ids_3olas <- a_bin_ao %>%
  filter(ola %in% c(1,3,6)) %>%
  count(idencuesta) %>%
  filter(n == 3) %>%
  pull(idencuesta)

dt_bal_ao <- a_bin_ao %>%
  filter(idencuesta %in% ids_3olas) %>%
  select(idencuesta, ola, all_of(items), n_member_ao)

cat("N personas:", length(ids_3olas), "| N obs:", nrow(dt_bal_ao), "\n")

cat("\nDistribución n_dominios (active_only) en panel balanceado:\n")
print(table(dt_bal_ao$n_member_ao))

# ------------------------------------------------------------------------------
# 4. CLASIFICACIÓN LMM (active_only)
# ------------------------------------------------------------------------------
cat("\n--- 4. CLASIFICACIÓN LMM (active_only: dt_states_cov_active_only.rds) ---\n")
dt_lmm_ao <- readRDS(here("data/dt_states_cov_active_only.rds")) %>%
  mutate(idencuesta = as.character(idencuesta))

cat("Tamaño de clases:\n")
print(dt_lmm_ao %>% count(position) %>% mutate(pct=round(100*n/sum(n),1)))

cat("\nn_member_ao MEDIO por clase (membresías activas observadas):\n")
dt_lmm_ao_bin <- dt_lmm_ao %>%
  left_join(
    dt_bal_ao %>% mutate(idencuesta=as.character(idencuesta)),
    by = c("idencuesta","ola")
  )
print(dt_lmm_ao_bin %>%
  group_by(position) %>%
  summarise(
    mean_n_ao   = round(mean(n_member_ao, na.rm=TRUE), 2),
    median_n    = median(n_member_ao, na.rm=TRUE),
    pct_0 = round(mean(n_member_ao==0, na.rm=TRUE)*100, 1),
    pct_1 = round(mean(n_member_ao==1, na.rm=TRUE)*100, 1),
    pct_2plus = round(mean(n_member_ao>=2, na.rm=TRUE)*100, 1)
  ))

cat("\nProporción observada c12==3 por clase y dominio:\n")
print(dt_lmm_ao_bin %>%
  group_by(position) %>%
  summarise(across(all_of(items), ~ round(mean(.x, na.rm=TRUE), 3))))

# ------------------------------------------------------------------------------
# 5. TRUST POR CLASE
# ------------------------------------------------------------------------------
cat("\n--- 5. TRUST POR CLASE (active_only) ---\n")
dt_anal <- readRDS(here("data/dt_analysis.rds")) %>%
  rename(idencuesta=id) %>%
  mutate(idencuesta=as.character(idencuesta)) %>%
  select(idencuesta, ola, trust, trust_nh)

dt_trust_ao <- dt_lmm_ao %>%
  left_join(dt_anal, by=c("idencuesta","ola"))
print(dt_trust_ao %>%
  group_by(position) %>%
  summarise(
    n = n(),
    trust_gen = round(mean(trust,    na.rm=TRUE)*100, 1),
    trust_nh  = round(mean(trust_nh, na.rm=TRUE)*100, 1)
  ))

# ------------------------------------------------------------------------------
# 6. COMPARACIÓN DIRECTA: ¿QUIÉNES CAMBIAN DE CLASE?
# ------------------------------------------------------------------------------
cat("\n--- 6. COMPARACIÓN any_member vs active_only ---\n")

dt_lmm_main <- readRDS(here("data/dt_states_cov.rds")) %>%
  mutate(idencuesta = as.character(idencuesta)) %>%
  select(idencuesta, ola, pos_main = position)

dt_cross <- inner_join(dt_lmm_main, dt_lmm_ao %>% select(idencuesta, ola, pos_ao=position),
                       by = c("idencuesta","ola"))

cat("\nTabla cruzada (filas=any_member, columnas=active_only):\n")
print(table(dt_cross$pos_main, dt_cross$pos_ao,
            dnn = c("any_member","active_only")))

cat("\n% flujos por fila:\n")
print(round(prop.table(table(dt_cross$pos_main, dt_cross$pos_ao), margin=1)*100, 1))

# --- Foco en el flujo más sospechoso: alpha(main) → gamma(ao) ---
cat("\n--- FOCO: alpha(any_member) → gamma(active_only) ---\n")
ids_ag <- dt_cross %>%
  filter(pos_main=="alpha", pos_ao=="gamma") %>%
  pull(idencuesta) %>% unique()

cat("N personas en este flujo:", length(ids_ag), "\n")

cat("\nSus c12 RAW (reales, sin recodificar) — media por dominio:\n")
a_raw_ag <- a %>%
  mutate(idencuesta = as.character(idencuesta)) %>%
  filter(idencuesta %in% ids_ag) %>%
  select(idencuesta, ola, all_of(items))

print(a_raw_ag %>%
  summarise(across(all_of(items), ~ {
    tbl <- table(.x)
    paste0("1=", tbl["1"], " 2=", tbl["2"], " 3=", tbl["3"])
  })) %>%
  pivot_longer(everything(), names_to="item", values_to="distribucion") %>%
  mutate(dominio = dominios[match(item, items)]) %>%
  select(dominio, distribucion))

cat("\nMean c12 (raw, 1-3) por dominio en este grupo:\n")
print(a_raw_ag %>%
  summarise(across(all_of(items), ~ round(mean(.x, na.rm=TRUE), 2))) %>%
  pivot_longer(everything(), names_to="item", values_to="mean_c12") %>%
  mutate(dominio = dominios[match(item, items)]) %>%
  arrange(desc(mean_c12)) %>%
  select(dominio, mean_c12))

cat("\nDistribución de n_member_ao y n_member_any para este grupo:\n")
a_raw_ag2 <- a_raw_ag %>%
  mutate(
    n_any = rowSums(across(all_of(items), ~ as.integer(!is.na(.x) & .x >= 2)), na.rm=TRUE),
    n_ao  = rowSums(across(all_of(items), ~ as.integer(!is.na(.x) & .x == 3)), na.rm=TRUE)
  )
cat("n_member_ao (active_only):\n")
print(table(a_raw_ag2$n_ao))
cat("n_member_any (any_member):\n")
print(table(a_raw_ag2$n_any))

cat("\n=== FIN diag_active_only.R ===\n")
