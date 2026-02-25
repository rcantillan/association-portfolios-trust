# ==============================================================================
# 04_precarity_transitions.R — "Structural precarity patterned by resources"
# ------------------------------------------------------------------------------
# Goal:
#   Show that exits from γ are patterned by resources (education/employment).
#
# Inputs:
#   data/dt_states.rds
# Outputs (output/):
#   - exit_gamma_models.txt
#   - exit_gamma_predicted.csv
#
# Design:
#   Create transitions between consecutive waves (t -> t+1).
#   Define exit_gamma = 1 if position_t == gamma and position_t1 != gamma.
#   Estimate logistic regression (or multinomial if you prefer) with resources.
# ==============================================================================

source(here::here("code", "00_setup.R"))

stop_if_missing(c(here::here("data", "dt_states.rds")))
dt <- readRDS(here::here("data", "dt_states.rds")) |> data.table::as.data.table()

dt[, position := factor(position, levels = c("alpha","beta","gamma"))]
data.table::setorder(dt, id, ola)

# Build transitions (t -> t+1) within id
dt[, position_next := data.table::shift(position, type = "lead"), by = id]
dt[, ola_next := data.table::shift(ola, type = "lead"), by = id]

trans <- dt[!is.na(position_next)]

# Exit from gamma
trans[, exit_gamma := as.integer(position == "gamma" & position_next != "gamma")]

# Restrict to risk set: those in gamma at t
risk <- trans[position == "gamma"]

# Covariates: education + employment (+ age/gender as baseline)
risk[, woman := as.integer(to01(woman))]
if (!is.numeric(risk$education)) risk[, education := factor(education)]
if (!is.numeric(risk$employment)) risk[, employment := factor(employment)]

m_exit <- glm(
  exit_gamma ~ education + employment + edad + woman + ola,
  data = risk,
  family = binomial(link = "logit")
)

# Predicted probabilities by education x employment (for a clean SSR figure/table)
# Build grid (if factors, keep all levels)
grid <- expand.grid(
  education = if (is.factor(risk$education)) levels(risk$education) else quantile(risk$education, probs = c(.25,.5,.75), na.rm = TRUE),
  employment = if (is.factor(risk$employment)) levels(risk$employment) else quantile(risk$employment, probs = c(.25,.5,.75), na.rm = TRUE),
  edad = median(risk$edad, na.rm = TRUE),
  woman = 0,
  ola = levels(factor(risk$ola))[1]
)

pred <- predict(m_exit, newdata = grid, type = "response")
out <- cbind(grid, pred_exit_gamma = pred)

readr::write_csv(as.data.frame(out), here::here("output", "exit_gamma_predicted.csv"))

sink(here::here("output", "exit_gamma_models.txt"))
cat("Logit model: exit from gamma (risk set: gamma at t)\n")
print(summary(m_exit))
sink()
