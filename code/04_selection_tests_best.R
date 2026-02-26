# ==============================================================================
# 04_selection_tests_best.R — Selection & temporal falsification (SSR block)
# ------------------------------------------------------------------------------
# Best version uses posterior uncertainty:
#  - Soft state: p_gamma (0..1)
#  - Soft entry: d_p_gamma = p_gamma(t+1) - p_gamma(t)
#  - Hard entry: enter_gamma_hard = I(position != gamma & next == gamma)
# Tests:
#  (i) baseline trust -> entry to gamma (soft + hard)
#  (ii) entry to gamma -> change in trust (Δtrust)
#  (iii) placebo lead: gamma_next predicts trust_t controlling gamma_now
# Outputs:
#  - selection_tests_table.tex (compact)
#  - CSV with key coefficients
# ==============================================================================

source(here::here("code","00_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(readr)
})

dt <- readRDS(here::here("data","dt_states.rds"))
dt <- as.data.table(dt)
setorder(dt, id, ola)

need <- c("id","ola","trust","trust_nh","position","p_gamma","edad","woman","education")
miss <- setdiff(need, names(dt))
if (length(miss) > 0) stop("Missing columns in dt_states.rds: ", paste(miss, collapse=", "))

# Leads
dt[, pos_next     := shift(position, -1), by=id]
dt[, trust_next   := shift(trust, -1),    by=id]
dt[, trust_nh_next:= shift(trust_nh, -1), by=id]
dt[, p_gamma_next := shift(p_gamma, -1),  by=id]

# Keep transition rows t -> t+1
tr <- dt[!is.na(pos_next)]

# Soft and hard indicators
tr[, gamma_now  := p_gamma]
tr[, gamma_next := p_gamma_next]
tr[, d_p_gamma  := p_gamma_next - p_gamma]  # soft entry (continuous)

tr[, enter_gamma_hard := as.integer(position != "gamma" & pos_next == "gamma")]
tr[, gamma_now_hard   := as.integer(position == "gamma")]
tr[, gamma_next_hard  := as.integer(pos_next  == "gamma")]

# Changes in trust (binary change -1/0/1)
tr[, d_trust    := trust_next - trust]
tr[, d_trust_nh := trust_nh_next - trust_nh]

# ---- Controls ----
# employment is often NA in your pipeline; we include only if it has variation.
ctrl <- c("edad","woman","education")
if ("employment" %in% names(tr) && any(!is.na(tr$employment)) && length(unique(na.omit(tr$employment))) > 1) {
  ctrl <- c(ctrl, "employment")
}

ctrl_str <- paste(ctrl, collapse = " + ")

# (i) Baseline trust -> entry to gamma
# (i-a) Soft: d_p_gamma ~ trust + gamma_now + controls + wave FE
f1_soft <- as.formula(paste0("d_p_gamma ~ trust + gamma_now + factor(ola) + ", ctrl_str))
m1_soft <- feols(f1_soft, data = tr, cluster = "id")

# (i-b) Hard: enter_gamma_hard ~ trust + gamma_now_hard + controls + wave FE (logit)
f1_hard <- as.formula(paste0("enter_gamma_hard ~ trust + gamma_now_hard + factor(ola) + ", ctrl_str))
m1_hard <- feglm(f1_hard, data = tr, family = "binomial", cluster = "id")

# (ii) Entry to gamma -> change in trust
# (ii-a) Soft: d_trust ~ d_p_gamma + trust + gamma_now + controls + wave FE
f2_soft <- as.formula(paste0("d_trust ~ d_p_gamma + trust + gamma_now + factor(ola) + ", ctrl_str))
m2_soft <- feols(f2_soft, data = tr, cluster = "id")

# (ii-b) Hard: d_trust ~ enter_gamma_hard + trust + gamma_now_hard + controls + wave FE
f2_hard <- as.formula(paste0("d_trust ~ enter_gamma_hard + trust + gamma_now_hard + factor(ola) + ", ctrl_str))
m2_hard <- feols(f2_hard, data = tr, cluster = "id")

# (iii) Placebo lead: gamma_next predicts trust_t controlling gamma_now (+ controls + wave FE)
# (iii-a) Soft placebo
f3_soft <- as.formula(paste0("trust ~ gamma_next + gamma_now + factor(ola) + ", ctrl_str))
m3_soft <- feglm(f3_soft, data = tr, family = "binomial", cluster = "id")

# (iii-b) Hard placebo
f3_hard <- as.formula(paste0("trust ~ gamma_next_hard + gamma_now_hard + factor(ola) + ", ctrl_str))
m3_hard <- feglm(f3_hard, data = tr, family = "binomial", cluster = "id")

# ---- Export table (compact SSR) ----
dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)

etable(
  "Baseline trust -> entry (soft)" = m1_soft,
  "Baseline trust -> entry (hard logit)" = m1_hard,
  "Entry -> Δtrust (soft)" = m2_soft,
  "Entry -> Δtrust (hard)" = m2_hard,
  "Placebo lead (soft logit)" = m3_soft,
  "Placebo lead (hard logit)" = m3_hard,
  file = here::here("output","selection_tests_table.tex"),
  tex = TRUE,
  replace = TRUE
)

# ---- Export key coefficients to CSV (for writing Results quickly) ----
grab <- function(model, term){
  co <- coef(model)
  se <- se(model)
  if (!(term %in% names(co))) return(data.table(term=term, estimate=NA_real_, se=NA_real_))
  data.table(term=term, estimate=unname(co[term]), se=unname(se[term]))
}

key <- rbindlist(list(
  cbind(test="m1_soft", grab(m1_soft, "trust")),
  cbind(test="m1_hard", grab(m1_hard, "trust")),
  cbind(test="m2_soft", grab(m2_soft, "d_p_gamma")),
  cbind(test="m2_hard", grab(m2_hard, "enter_gamma_hard")),
  cbind(test="m3_soft", grab(m3_soft, "gamma_next")),
  cbind(test="m3_hard", grab(m3_hard, "gamma_next_hard"))
), fill=TRUE)

write_csv(key, here::here("output","selection_tests_keycoef.csv"))

saveRDS(list(m1_soft=m1_soft, m1_hard=m1_hard, m2_soft=m2_soft, m2_hard=m2_hard, m3_soft=m3_soft, m3_hard=m3_hard),
        here::here("output","selection_tests_models.rds"))

message("DONE: output/selection_tests_table.tex + selection_tests_keycoef.csv")