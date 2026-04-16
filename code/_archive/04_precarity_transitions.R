# ==============================================================================
# 04_precarity_transitions.R — Structural precarity: resources -> exits from γ
# ------------------------------------------------------------------------------
# Input:
#   data/dt_states.rds
# Outputs:
#   output/exit_gamma_logit_summary.txt
#   output/exit_gamma_predicted.csv
# ==============================================================================

source(here::here("code", "00_setup.R"))

stop_if_missing(c(here::here("data", "dt_states.rds")))
dt <- readRDS(here::here("data", "dt_states.rds")) |> as_tibble()

dt <- dt %>%
  arrange(id, ola) %>%
  mutate(position = factor(position, levels = c("alpha","beta","gamma")))

# Build t -> t+1 transitions (balanced panel already)
trans <- dt %>%
  group_by(id) %>%
  mutate(
    position_next = lead(position),
    ola_next = lead(ola)
  ) %>%
  ungroup() %>%
  filter(!is.na(position_next))

# Risk set: those in γ at time t
risk_gamma <- trans %>%
  filter(position == "gamma") %>%
  mutate(exit_gamma = as.integer(position_next != "gamma"))

# Minimal resource predictors (education/employment) + controls
# (Adjust education recode later if you want categories)
m_exit <- glm(exit_gamma ~ education + employment + edad + woman + factor(ola),
              data = risk_gamma, family = binomial(), na.action = na.omit)

# Save model summary
out_txt <- capture.output(summary(m_exit))
writeLines(out_txt, con = here::here("output", "exit_gamma_logit_summary.txt"))

# Predicted exit probabilities by education/employment (simple grid)
grid <- expand.grid(
  education  = sort(unique(risk_gamma$education[!is.na(risk_gamma$education)])),
  employment = sort(unique(risk_gamma$employment[!is.na(risk_gamma$employment)])),
  edad = median(risk_gamma$edad, na.rm = TRUE),
  woman = 1,
  ola = unique(risk_gamma$ola)[1]
)

grid$pred_exit <- predict(m_exit, newdata = grid, type = "response")
readr::write_csv(as_tibble(grid), here::here("output", "exit_gamma_predicted.csv"))

print(m_exit)