# ==============================================================================
# any_member / 02_lmm.R
# Latent Markov Model K=3 sobre membresías (codificación ANY_MEMBER)
# Input:  output/dt_panel.rds
# Output: output/lmm_fit.csv
#         output/lmm_profiles.csv
#         output/lmm_transitions.csv
#         output/dt_states.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(LMest); library(data.table); library(here)
})

ROOT    <- here::here()
DIR_OUT <- file.path(ROOT, "any_member", "output")

ITEMS   <- c(paste0("c12_0", 1:8), "c12_09")
DOMAINS <- c("nhg","religious","political","union",
              "professional","charity","sport","student","other")
K_GRID  <- 1:5
K_USE   <- 3L
SEED    <- 123L

cat("=== 02_lmm.R | any_member | K=", K_USE, "===\n\n")

# ------------------------------------------------------------------------------
# 1. Cargar datos
# ------------------------------------------------------------------------------
dt <- readRDS(file.path(DIR_OUT, "dt_panel.rds"))
cat("Datos cargados:", nrow(dt), "obs |", n_distinct(dt$idencuesta), "personas\n\n")

# Variables LMM: ids, olas, items + covariables de transición
# Covariables: woman, edad categórica, educación
dt <- dt %>%
  mutate(
    mujer       = woman,
    edad_cat    = cut(age,
                      breaks  = c(17,24,34,44,54,64,120),
                      labels  = c("18_24","25_34","35_44","45_54","55_64","65+"),
                      right   = TRUE),
    nivel_educ  = case_when(
      education == 1L ~ "basica",
      education == 2L ~ "media",
      education == 3L ~ "tecnica",
      education == 4L ~ "univers",
      TRUE            ~ NA_character_
    )
  ) %>%
  filter(!is.na(mujer), !is.na(edad_cat), !is.na(nivel_educ))

# Re-balancear tras drop de NAs en covariables
ids_ok <- dt %>% count(idencuesta) %>% filter(n == 3L) %>% pull(idencuesta)
dt <- filter(dt, idencuesta %in% ids_ok)
cat("Tras drop NA (covariables LMM):", n_distinct(dt$idencuesta), "personas |",
    nrow(dt), "obs\n\n")

# ------------------------------------------------------------------------------
# 2. Selección K con lmestSearch
# ------------------------------------------------------------------------------
set.seed(SEED)
cat("Corriendo lmestSearch K=1 a", max(K_GRID), "...\n")
mod_sel <- lmestSearch(
  responsesFormula = c12_01 + c12_02 + c12_03 + c12_04 +
                     c12_05 + c12_06 + c12_07 + c12_08 + c12_09 ~ NULL,
  latentFormula    = ~ mujer + edad_cat + nivel_educ,
  index            = c("idencuesta", "ola"),
  data             = dt,
  output           = TRUE,
  out_se           = FALSE,
  version          = "categorical",
  paramLatent      = "multilogit",
  k                = K_GRID,
  seed             = SEED
)

fit_tbl <- tibble(
  K      = K_GRID,
  logLik = as.numeric(mod_sel$lkv),
  AIC    = as.numeric(mod_sel$Aic),
  BIC    = as.numeric(mod_sel$Bic)
) %>% mutate(delta_BIC = BIC - lag(BIC))

cat("\n=== FIT TABLE ===\n")
print(fit_tbl)
write_csv(fit_tbl, file.path(DIR_OUT, "lmm_fit.csv"))

K_star <- fit_tbl$K[which.min(fit_tbl$BIC)]
cat("\nBIC selecciona K =", K_star, "| Usamos K =", K_USE, "\n\n")

# ------------------------------------------------------------------------------
# 3. Extraer modelo K = K_USE
# ------------------------------------------------------------------------------
mod3 <- mod_sel$out.single[[K_USE]]

Psi <- mod3$Psi   # item x state x category
Pi  <- mod3$Pi    # transition matrix
V   <- mod3$V     # posteriors

ids   <- sort(unique(dt$idencuesta))
waves <- sort(unique(dt$ola))
n     <- length(ids)
TT    <- length(waves)
r     <- length(ITEMS)

# Reshape V → n x K x TT
reshape_V <- function(V, n, k, TT) {
  perms <- list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
  best <- NULL; best_err <- Inf
  for (p in perms) {
    Vp <- aperm(V, p)
    if (!all(dim(Vp) == c(n,k,TT))) next
    err <- max(abs(apply(Vp, c(1,3), sum) - 1), na.rm = TRUE)
    if (err < best_err) { best_err <- err; best <- Vp }
  }
  if (is.null(best)) stop("No se pudo reshape V a (n,k,TT)")
  best
}
V_std <- reshape_V(V, n, K_USE, TT)

# P(Y=1|state): r x K
Psi_std <- NULL
for (p in list(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))) {
  Pp <- aperm(Psi, p)
  if (all(dim(Pp) == c(r, K_USE, 2))) { Psi_std <- Pp; break }
}
prob1 <- Psi_std[,,2]  # r x K
rownames(prob1) <- ITEMS

# ------------------------------------------------------------------------------
# 4. Mapear estados → α / β / γ
# ------------------------------------------------------------------------------
exp_count   <- colSums(prob1, na.rm = TRUE)
state_entr  <- sapply(1:K_USE, function(s) {
  p <- prob1[,s] / (sum(prob1[,s]) + 1e-12)
  -sum(p * log(p + 1e-12)) / log(r)
})
names(exp_count) <- names(state_entr) <- paste0("State_", 1:K_USE)

alpha_s <- which.min(exp_count)
rem     <- setdiff(1:K_USE, alpha_s)
gamma_s <- rem[which.max(state_entr[rem])]
beta_s  <- setdiff(1:K_USE, c(alpha_s, gamma_s))

cat("=== STATE MAPPING ===\n")
cat(sprintf("  α (Isolate)  = State %d | exp_count=%.3f | entropy=%.3f\n",
            alpha_s, exp_count[alpha_s], state_entr[alpha_s]))
cat(sprintf("  β (Closed)   = State %d | exp_count=%.3f | entropy=%.3f\n",
            beta_s,  exp_count[beta_s],  state_entr[beta_s]))
cat(sprintf("  γ (Bridging) = State %d | exp_count=%.3f | entropy=%.3f\n\n",
            gamma_s, exp_count[gamma_s], state_entr[gamma_s]))

# Verificación sustantiva
n_high_g <- sum(prob1[, gamma_s] >= 0.30)
n_high_b <- sum(prob1[, beta_s]  >= 0.30)
cat(sprintf("  Dominios con Pr>=0.30: α=%d | β=%d | γ=%d (esperado γ>=4)\n\n",
            sum(prob1[,alpha_s]>=0.30), n_high_b, n_high_g))
if (n_high_g < 3)
  warning("γ tiene solo ", n_high_g, " dominios con Pr>=0.30 — verificar mapeo")

# ------------------------------------------------------------------------------
# 5. Perfiles y transiciones
# ------------------------------------------------------------------------------
prof_df <- as.data.frame(round(prob1, 3))
colnames(prof_df) <- paste0("State_", 1:K_USE)
prof_df$item   <- ITEMS
prof_df$domain <- DOMAINS
prof_df$state_alpha <- prof_df[[paste0("State_", alpha_s)]]
prof_df$state_beta  <- prof_df[[paste0("State_", beta_s)]]
prof_df$state_gamma <- prof_df[[paste0("State_", gamma_s)]]

cat("=== PERFILES (Pr membresía | clase) ===\n")
print(prof_df %>% select(domain, state_alpha, state_beta, state_gamma))

write_csv(prof_df %>% select(domain, state_alpha, state_beta, state_gamma),
          file.path(DIR_OUT, "lmm_profiles.csv"))

# Transición
Pi_avg <- if (length(dim(Pi))==3) apply(Pi, c(1,2), mean) else Pi
ord    <- c(alpha_s, beta_s, gamma_s)
Pi_r   <- Pi_avg[ord, ord, drop=FALSE]
rownames(Pi_r) <- colnames(Pi_r) <- c("alpha","beta","gamma")

cat("\n=== MATRIZ DE TRANSICIÓN (media de olas) ===\n")
print(round(Pi_r, 3))
write_csv(as.data.frame(Pi_r) %>% rownames_to_column("from"),
          file.path(DIR_OUT, "lmm_transitions.csv"))

# ------------------------------------------------------------------------------
# 6. Asignaciones modales + posteriors
# ------------------------------------------------------------------------------
label_s <- function(s) dplyr::case_when(
  s == alpha_s ~ "alpha", s == gamma_s ~ "gamma", TRUE ~ "beta"
)

post_long <- rbindlist(lapply(seq_len(TT), function(t) {
  pp <- V_std[,,t]
  data.table(
    idencuesta = ids, ola = waves[t],
    p_alpha = pp[, alpha_s],
    p_beta  = pp[, beta_s],
    p_gamma = pp[, gamma_s],
    p_max   = apply(pp, 1, max)
  )
}))

modal_dt <- rbindlist(lapply(seq_len(TT), function(t) {
  Uhat <- max.col(V_std[,,t], ties.method="first")
  data.table(idencuesta=ids, ola=waves[t],
             position=label_s(Uhat))
}))

dt_states <- merge(
  merge(as.data.table(dt), post_long, by=c("idencuesta","ola"), all.x=TRUE),
  modal_dt,  by=c("idencuesta","ola"), all.x=TRUE
)
dt_states[, position := factor(position, levels=c("alpha","beta","gamma"))]

cat("\n=== TAMAÑO DE CLASES ===\n")
print(dt_states[, .N, by=position][order(position)][
  , pct := round(100*N/sum(N),1)])

cat("\n=== CALIDAD DE CLASIFICACIÓN ===\n")
print(dt_states[, .(mean_pmax=round(mean(p_max),3),
                     share_60  =round(mean(p_max>=0.60),3),
                     n=.N)])

saveRDS(as_tibble(dt_states), file.path(DIR_OUT, "dt_states.rds"))
cat("\n✓ Guardado:", file.path(DIR_OUT, "dt_states.rds"), "\n")
cat("\n[02_lmm.R DONE]\n")
