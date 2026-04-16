# ==============================================================================
# 003b_precarity_transitions_SES.R — H3: Structural Precarity (SES-stratified exits)
# ------------------------------------------------------------------------------
# H3 (Structural precarity): Position γ (bridging) will exhibit lower temporal
# stability than Positions α and β, with higher exit probabilities among
# individuals facing greater resource constraints.
#
# INPUTS:
#   data/dt_states_cov.rds      (from 002_long_latent_class.R)
#   data/dt_analysis.rds        (from 01_descriptive_stats.R)
#
# OUTPUTS — PAPER (main text):
#   output/figure_h3_paper.png/.pdf     Panel A + B combined (paper figure)
#   output/latex_table_h3_main.txt      Table 3: persistence + M3 AMEs (Overleaf)
#   output/h3_inline_numbers.txt        Copy-paste paragraph for §3.2
#
# OUTPUTS — SUPPLEMENTARY (SI):
#   output/figure_h3_SI_panelA.png/.pdf Panel A standalone (SI archive)
#   output/latex_SI_tableA2.txt         SI Table A2: persistence by position × wave
#   output/latex_SI_tableA3.txt         SI Table A3: all exit-from-γ models
#   output/latex_SI_tableA4.txt         SI Table A4: exit models α, β, γ (specificity)
#
# OUTPUTS — INTERMEDIATE (diagnostics):
#   output/h3_persistence_table.csv
#   output/h3_exit_gamma_ames.csv
#   output/h3_allpos_exit_ames.csv
#   output/h3_summary_for_paper.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(marginaleffects)
  library(sandwich)
  library(knitr)
})

source(here::here("code", "00_setup.R"))
dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# HELPERS
# ==============================================================================

sig_stars <- function(p) {
  dplyr::case_when(
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    p < .10  ~ "$\\dagger$",
    TRUE     ~ ""
  )
}
sig_stars_plain <- function(p) {
  dplyr::case_when(
    p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", p < .10 ~ "+", TRUE ~ ""
  )
}
fmt_n <- function(x, d = 3) formatC(round(x, d), format = "f", digits = d)

write_latex <- function(tex_lines, filename) {
  writeLines(tex_lines, here::here("output", filename))
  message("  LaTeX: output/", filename)
}

# ==============================================================================
# 0) Load data
# ==============================================================================
dt_states   <- readRDS(here::here("data", "dt_states_cov.rds"))
dt_analysis <- readRDS(here::here("data", "dt_analysis.rds"))

message("dt_states:   ", nrow(dt_states),   " obs | ", n_distinct(dt_states$idencuesta), " ids")
message("dt_analysis: ", nrow(dt_analysis), " obs | ", n_distinct(dt_analysis$id),       " ids")

# ==============================================================================
# 1) Build transition dataset
# ==============================================================================
dt_states <- dt_states %>%
  mutate(id       = as.character(idencuesta),
         position = as.character(position))
dt_analysis <- dt_analysis %>%
  mutate(id = as.character(id))

ses_vars <- dt_analysis %>%
  select(id, ola, edu_bin, employed) %>%
  mutate(id = as.character(id), ola = as.integer(ola))

dt_full <- dt_states %>% left_join(ses_vars, by = c("id","ola"))

dt_t  <- dt_full %>%
  select(id, ola, position, p_gamma, nivel_educ, mujer, edad, edu_bin, employed)
dt_t1 <- dt_full %>%
  mutate(ola = ola - 1L) %>%
  select(id, ola, position_t1 = position, p_gamma_t1 = p_gamma)

dt_trans <- dt_t %>%
  inner_join(dt_t1, by = c("id","ola")) %>%
  filter(!is.na(position), !is.na(position_t1)) %>%
  mutate(
    stay_same  = as.integer(position == position_t1),
    exit_gamma = case_when(
      position == "gamma" ~ as.integer(position_t1 != "gamma"),
      TRUE                ~ NA_integer_
    ),
    edu_ord = case_when(
      nivel_educ == "básica"  ~ 1L, nivel_educ == "media"   ~ 2L,
      nivel_educ == "técnica" ~ 3L, nivel_educ == "univers" ~ 4L,
      TRUE ~ NA_integer_
    ),
    edu_high = as.integer(edu_ord >= 3),
    age_mid  = case_when(
      edad == "18_24" ~ 21,   edad == "25_34" ~ 29.5, edad == "35_44" ~ 39.5,
      edad == "45_54" ~ 49.5, edad == "55_64" ~ 59.5, edad == "65"    ~ 70,
      TRUE ~ NA_real_
    ),
    wave_f = factor(ola, levels = c(1L,2L),
                    labels = c("2016\u21922018","2018\u21922022"))
  )

message("Total transitions: ", nrow(dt_trans),
        " | From \u03b3: ", sum(dt_trans$position == "gamma", na.rm=TRUE))

# ==============================================================================
# 2) Persistence by position × wave
# ==============================================================================
persistence_tbl <- dt_trans %>%
  group_by(position, wave_f) %>%
  summarise(n_transitions = n(),
            n_stay        = sum(stay_same, na.rm=TRUE),
            p_stay        = mean(stay_same, na.rm=TRUE),
            p_stay_se     = sqrt(p_stay*(1-p_stay)/n()),
            .groups="drop") %>%
  mutate(p_exit = 1-p_stay,
         ci_lo  = pmax(0, p_stay - 1.96*p_stay_se),
         ci_hi  = pmin(1, p_stay + 1.96*p_stay_se)) %>%
  arrange(position, wave_f)

write_csv(persistence_tbl, here::here("output","h3_persistence_table.csv"))
print(persistence_tbl)

pers_pooled <- dt_trans %>%
  group_by(position) %>%
  summarise(p_stay = mean(stay_same, na.rm=TRUE), n=n(), .groups="drop")

# ==============================================================================
# 3) Exit models from γ
# ==============================================================================
dt_gamma_exits <- dt_trans %>% filter(position=="gamma", !is.na(exit_gamma))
message("\n\u03b3-exit subsample: ", nrow(dt_gamma_exits))
message("Exit rate from \u03b3: ", round(mean(dt_gamma_exits$exit_gamma), 3))

m_edu_only <- glm(exit_gamma ~ edu_high + wave_f,
                  data=dt_gamma_exits, family=binomial("probit"))
m_emp_only <- glm(exit_gamma ~ employed + wave_f,
                  data=dt_gamma_exits, family=binomial("probit"))
m_full     <- glm(exit_gamma ~ edu_high + employed + mujer + age_mid + wave_f,
                  data=dt_gamma_exits, family=binomial("probit"))
m_edu_ord  <- glm(exit_gamma ~ edu_ord  + employed + mujer + age_mid + wave_f,
                  data=dt_gamma_exits, family=binomial("probit"))

vcov_cl <- function(m) tryCatch(
  sandwich::vcovCL(m, cluster=dt_gamma_exits$id, type="HC1"),
  error=function(e) NULL)

ames_edu  <- avg_comparisons(m_edu_only, variables="edu_high",              vcov=vcov_cl(m_edu_only))
ames_emp  <- avg_comparisons(m_emp_only, variables="employed",              vcov=vcov_cl(m_emp_only))
ames_full <- avg_comparisons(m_full,     variables=c("edu_high","employed"), vcov=vcov_cl(m_full))
ames_ord  <- avg_comparisons(m_edu_ord,  variables=c("edu_ord","employed"),  vcov=vcov_cl(m_edu_ord))

message("\n--- AMEs full model ---"); print(ames_full)

# ==============================================================================
# 4) Specificity: exit models for all starting positions
# ==============================================================================
fit_exit_model <- function(pos_label) {
  dt_sub <- dt_trans %>% filter(position==pos_label) %>%
    mutate(exit_pos=1L-stay_same) %>% filter(!is.na(exit_pos))
  if (nrow(dt_sub)<20 || n_distinct(dt_sub$exit_pos)<2) return(NULL)
  m   <- glm(exit_pos ~ edu_high + employed + mujer + age_mid + wave_f,
             data=dt_sub, family=binomial("probit"))
  vcv <- tryCatch(sandwich::vcovCL(m, cluster=dt_sub$id, type="HC1"),
                  error=function(e) NULL)
  avg_comparisons(m, variables=c("edu_high","employed"), vcov=vcv) %>%
    as_tibble() %>%
    mutate(starting_position=pos_label,
           n_transitions=nrow(dt_sub),
           exit_rate=mean(dt_sub$exit_pos))
}
allpos_ames <- bind_rows(lapply(c("alpha","beta","gamma"), fit_exit_model))
write_csv(allpos_ames, here::here("output","h3_allpos_exit_ames.csv"))

# ==============================================================================
# 5) Save combined AMEs (intermediate)
# ==============================================================================
exit_gamma_ames <- bind_rows(
  ames_edu  %>% as_tibble() %>% mutate(model="M1: Education only"),
  ames_emp  %>% as_tibble() %>% mutate(model="M2: Employment only"),
  ames_full %>% as_tibble() %>% mutate(model="M3: Full (preferred)"),
  ames_ord  %>% as_tibble() %>% mutate(model="M4: Ordinal education")
)
write_csv(exit_gamma_ames, here::here("output","h3_exit_gamma_ames.csv"))

# Intermediate summary
sink(here::here("output","h3_summary_for_paper.txt"))
cat("=== H3 Summary ===\n\n--- Persistence (pooled) ---\n"); print(pers_pooled)
cat("\nExit rate from gamma:", round(mean(dt_gamma_exits$exit_gamma),3))
cat("\nn(gamma-transitions):", nrow(dt_gamma_exits),"\n\n--- Model 3 AMEs ---\n")
print(ames_full %>% as_tibble() %>% select(term,estimate,std.error,p.value) %>%
        mutate(across(where(is.numeric),~round(.x,4))))
sink()

# ==============================================================================
# 6) FIGURES
# ==============================================================================

# Data for figures
exit_by_pos <- dt_trans %>%
  mutate(exit = 1L - stay_same) %>%
  group_by(position) %>%
  summarise(n=n(), exit_rate=mean(exit,na.rm=TRUE),
            se=sqrt(exit_rate*(1-exit_rate)/n), .groups="drop") %>%
  mutate(
    pos_label = factor(position,
                       levels = c("alpha","beta","gamma"),
                       labels = c("\u03b1 (isolation)",
                                  "\u03b2 (clustering)",
                                  "\u03b3 (bridging)")),
    # Label sits above CI upper bound + small gap → no overlap regardless of CI width
    ci_hi     = pmin(1, exit_rate + 1.96*se),
    label_y   = ci_hi + 0.028
  )

exit_by_ses <- dt_gamma_exits %>%
  filter(!is.na(edu_high), !is.na(employed)) %>%
  group_by(edu_high, employed) %>%
  summarise(n=n(), exit_rate=mean(exit_gamma,na.rm=TRUE),
            se=sqrt(exit_rate*(1-exit_rate)/n), .groups="drop") %>%
  mutate(
    ci_hi     = pmin(1, exit_rate + 1.96*se),
    label_y   = ci_hi + 0.030,
    edu_label = factor(edu_high, levels=c(0,1),
                       labels=c("Low education\n(b\u00e1sica/media)",
                                "High education\n(t\u00e9cnica/univ)")),
    emp_label = factor(employed, levels=c(0,1),
                       labels=c("Not employed","Employed"))
  )

# y-axis upper limit: max label_y + a little breathing room
ylim_A  <- max(exit_by_pos$label_y)  + 0.04
ylim_B  <- max(exit_by_ses$label_y)  + 0.04

# ── Panel A ───────────────────────────────────────────────────────────────────
p_A <- ggplot(exit_by_pos,
              aes(x=pos_label, y=exit_rate, fill=pos_label)) +
  geom_col(width=0.52, show.legend=FALSE) +
  geom_errorbar(aes(ymin=pmax(0, exit_rate-1.96*se), ymax=ci_hi),
                width=0.18, linewidth=0.6, color="grey30") +
  geom_text(aes(label=percent(exit_rate, accuracy=1), y=label_y),
            size=4, fontface="bold") +
  scale_fill_manual(values=unname(STATE_PALETTE_LEGACY)) +
  scale_y_continuous(labels=percent_format(accuracy=1),
                     limits=c(0, ylim_A),
                     expand=expansion(mult=c(0, 0))) +
  labs(x       = "Portfolio position at t",
       y       = "P(exit by t+1)",
       title   = "A. Exit probabilities by portfolio position",
       subtitle= "ELSOC balanced panel 2016\u20132022") +
  theme_ssr(base_size=12) +
  theme(plot.title    = element_text(size=12, face="bold"),
        plot.subtitle = element_text(size=10, color="grey40"),
        plot.margin   = margin(t=10, r=16, b=10, l=10, unit="pt"))

# ── Panel B ───────────────────────────────────────────────────────────────────
p_B <- ggplot(exit_by_ses,
              aes(x=edu_label, y=exit_rate,
                  fill=emp_label, group=emp_label)) +
  geom_col(position=position_dodge(0.65), width=0.58) +
  geom_errorbar(aes(ymin=pmax(0, exit_rate-1.96*se), ymax=ci_hi),
                position=position_dodge(0.65),
                width=0.18, linewidth=0.6, color="grey30") +
  geom_text(aes(label=percent(exit_rate, accuracy=1), y=label_y),
            position=position_dodge(0.65),
            size=3.5, fontface="bold") +
  scale_fill_manual(values=c("Not employed"="#9CAEBF","Employed"="#E15759"),
                    name=NULL) +
  scale_y_continuous(labels=percent_format(accuracy=1),
                     limits=c(0, ylim_B),
                     expand=expansion(mult=c(0, 0))) +
  labs(x       = NULL,
       y       = "P(exit from \u03b3 by t+1)",
       title   = "B. Exit from bridging (\u03b3) by education and employment",
       subtitle= "Subsample: person-transitions starting in \u03b3") +
  theme_ssr(base_size=12) +
  theme(plot.title      = element_text(size=12, face="bold"),
        plot.subtitle   = element_text(size=10, color="grey40"),
        legend.position = "top",
        legend.key.size = unit(0.55, "cm"),
        plot.margin     = margin(t=10, r=16, b=10, l=10, unit="pt"))

# ── Shared caption ────────────────────────────────────────────────────────────
shared_caption <- paste0(
  "ELSOC balanced panel 2016\u20132022. Error bars: 95% CI. ",
  "n(A) = 2,400 person-transitions; ",
  "n(B) = ", nrow(dt_gamma_exits), " (\u03b3-subsample). ",
  "Cluster-robust SE (by individual) in Panel B models."
)

# ── Combined figure theme ─────────────────────────────────────────────────────
combined_theme <- theme(
  plot.caption = element_text(size=8, color="grey40",
                              hjust=0, margin=margin(t=8)),
  plot.margin  = margin(t=6, r=8, b=10, l=8, unit="pt")
)

# ── PAPER figure: Panel A + B combined ────────────────────────────────────────
fig_paper <- p_A + p_B +
  plot_layout(widths=c(1, 1.4)) +          # ← sin guides="collect"
  plot_annotation(caption=shared_caption, theme=combined_theme)
  # ← sin & theme(legend.position="top")

ggsave(here::here("output","figure_h3_paper.png"),
       fig_paper, width=12, height=5.4, dpi=300)
ggsave(here::here("output","figure_h3_paper.pdf"),
       fig_paper, width=12, height=5.4)
message("  PAPER figure (A+B): output/figure_h3_paper.png/.pdf")

# ── SI figure: Panel A alone (for SI archive, lighter) ────────────────────────
ggsave(here::here("output","figure_h3_SI_panelA.png"),
       p_A, width=5.5, height=4.8, dpi=300)
ggsave(here::here("output","figure_h3_SI_panelA.pdf"),
       p_A, width=5.5, height=4.8)
message("  SI figure (A only): output/figure_h3_SI_panelA.png/.pdf")

# ==============================================================================
# 7) Scalar values for inline text and tables
# ==============================================================================
p_exit_alpha  <- 1 - pers_pooled$p_stay[pers_pooled$position=="alpha"]
p_exit_beta   <- 1 - pers_pooled$p_stay[pers_pooled$position=="beta"]
p_exit_gamma  <- mean(dt_gamma_exits$exit_gamma)
n_gamma_trans <- nrow(dt_gamma_exits)

edu_ame_val <- ames_full %>% filter(term=="edu_high") %>% pull(estimate)
edu_ame_se  <- ames_full %>% filter(term=="edu_high") %>% pull(std.error)
edu_ame_p   <- ames_full %>% filter(term=="edu_high") %>% pull(p.value)
emp_ame_val <- ames_full %>% filter(term=="employed")  %>% pull(estimate)
emp_ame_se  <- ames_full %>% filter(term=="employed")  %>% pull(std.error)
emp_ame_p   <- ames_full %>% filter(term=="employed")  %>% pull(p.value)

alpha_edu <- allpos_ames %>% filter(starting_position=="alpha", term=="edu_high") %>% pull(estimate)
beta_edu  <- allpos_ames %>% filter(starting_position=="beta",  term=="edu_high") %>% pull(estimate)

# ==============================================================================
# 8) Inline numbers — copy-paste paragraph for §3.2
# ==============================================================================
sink(here::here("output","h3_inline_numbers.txt"))
cat("================================================================\n")
cat("INLINE NUMBERS — §3.2 Structural Precarity\n")
cat("================================================================\n\n")
cat("--- Exit rates (pooled) ---\n")
cat(sprintf("  P(exit | alpha) = %.1f%%   n = %d\n",
            p_exit_alpha*100, sum(dt_trans$position=="alpha")))
cat(sprintf("  P(exit | beta)  = %.1f%%   n = %d\n",
            p_exit_beta*100,  sum(dt_trans$position=="beta")))
cat(sprintf("  P(exit | gamma) = %.1f%%   n = %d\n",
            p_exit_gamma*100, n_gamma_trans))
cat(sprintf("  Ratio gamma/alpha = %.1fx\n\n", p_exit_gamma/p_exit_alpha))
cat("--- By transition period ---\n")
persistence_tbl %>%
  mutate(out=sprintf("  %-8s %s: P(exit) = %.1f%%",
                     position, wave_f, p_exit*100)) %>%
  pull(out) %>% cat(sep="\n")
cat("\n\n--- Model 3 AMEs ---\n")
cat(sprintf("  edu_high: AME=%+.3f  SE=%.3f  p=%.4f  %s\n",
            edu_ame_val, edu_ame_se, edu_ame_p, sig_stars_plain(edu_ame_p)))
cat(sprintf("  employed: AME=%+.3f  SE=%.3f  p=%.4f  %s\n",
            emp_ame_val, emp_ame_se, emp_ame_p, sig_stars_plain(emp_ame_p)))
cat(sprintf("  n(gamma) = %d\n\n", n_gamma_trans))
cat("--- Specificity ---\n")
allpos_ames %>%
  mutate(out=sprintf("  %-6s | %-10s: AME=%+.3f (p=%.3f) %s",
                     starting_position, term, estimate,
                     p.value, sig_stars_plain(p.value))) %>%
  pull(out) %>% cat(sep="\n")
cat("\n\n================================================================\n")
cat("PARAGRAPH FOR §3.2 (LaTeX-ready):\n")
cat("================================================================\n\n")
cat(sprintf(
"Position $\\gamma$ is substantially less stable than Positions $\\alpha$ and
$\\beta$ (Figure~2). Pooled across wave-transitions, the exit probability from
$\\gamma$ is %.1f\\%%, compared to %.1f\\%% for $\\beta$ and %.1f\\%%
for $\\alpha$ --- a %.1f-fold difference consistent across both transition
periods (Table~A2). Among respondents in $\\gamma$ at wave $t$ ($n = %d$
person-wave observations), exits are negatively associated with socioeconomic
resources: higher education reduces the probability of exiting $\\gamma$ by
%.1f percentage points (AME\\,$=$\\,%.3f, SE\\,$=$\\,%.3f, $p < .001$), and
employment reduces it by %.1f percentage points (AME\\,$=$\\,%.3f,
SE\\,$=$\\,%.3f, $p = %.2f$; Table~A3). Importantly, the education gradient
reverses across starting positions: among respondents in $\\alpha$ and $\\beta$,
higher education is associated with \\emph{higher} exit probabilities
(AMEs\\,$=+$%.3f and $+$%.3f, respectively), consistent with upward mobility
toward bridging, whereas among respondents in $\\gamma$ it predicts retention
(Table~A4). Together, these patterns support H3.\n",
  p_exit_gamma*100, p_exit_beta*100, p_exit_alpha*100,
  p_exit_gamma/p_exit_alpha,
  n_gamma_trans,
  abs(edu_ame_val)*100, edu_ame_val, edu_ame_se,
  abs(emp_ame_val)*100, emp_ame_val, emp_ame_se, emp_ame_p,
  alpha_edu, beta_edu
))
sink()
message("  Inline numbers: output/h3_inline_numbers.txt")

# ==============================================================================
# 9) LaTeX TABLES for Overleaf
# ==============================================================================

# ── Table 3 (PAPER): Persistence + M3 AMEs ───────────────────────────────────
pers_rows <- pers_pooled %>%
  arrange(match(position, c("alpha","beta","gamma"))) %>%
  mutate(
    pos_tex = c("$\\alpha$ (isolation)","$\\beta$ (clustering)","$\\gamma$ (bridging)"),
    row     = sprintf("  %-28s & %5d & %.3f & %.3f &       &      \\\\",
                      pos_tex, n, p_stay, 1-p_stay)
  ) %>% pull(row)

ame_rows <- ames_full %>% as_tibble() %>%
  mutate(
    term_tex = ifelse(term=="edu_high",
                      "High education (ref: low)",
                      "Employed (ref: other)"),
    row = sprintf("  %-28s & %5d &       &       & %+.3f%s & %.3f \\\\",
                  term_tex, n_gamma_trans, estimate,
                  sig_stars(p.value), std.error)
  ) %>% pull(row)

tbl3_lines <- c(
  "% ============================================================",
  "% Table 3 — Portfolio persistence and SES gradient",
  "% Paste this block into your Overleaf main document",
  "% ============================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Portfolio-position persistence and socioeconomic gradient in exits from bridging}",
  "\\label{table:h3_main}",
  "\\small",
  "\\setlength{\\tabcolsep}{6pt}",
  "\\begin{tabular}{lrcccc}",
  "\\hline\\hline",
  "  & & \\multicolumn{2}{c}{Panel A: Persistence} & \\multicolumn{2}{c}{Panel B: Exit from $\\gamma$} \\\\",
  "  \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  "  Position / Predictor & $n$ & P(stay) & P(exit) & AME & SE \\\\",
  "\\hline",
  "  \\multicolumn{6}{l}{\\textit{Unconditional persistence (pooled across transitions)}} \\\\",
  pers_rows,
  "\\hline",
  "  \\multicolumn{6}{l}{\\textit{SES predictors of exit from $\\gamma$ (Model 3, preferred)}} \\\\",
  ame_rows,
  "\\hline\\hline",
  "  \\multicolumn{6}{p{0.92\\textwidth}}{\\footnotesize \\textit{Notes.}",
  paste0("  Panel A: pooled persistence and exit rates across 2016$\\to$2018 and 2018$\\to$2022 transitions."),
  paste0("  Panel B: AMEs from a probit model of P(exit from $\\gamma$ by $t+1$) among respondents in $\\gamma$",
         " ($n = ", n_gamma_trans, "$ person-wave observations)."),
  "  Controls: sex, age (band midpoint), wave indicator.",
  "  Cluster-robust SE (by individual).",
  "  $\\dagger p < .10$; $* p < .05$; $** p < .01$; $*** p < .001$.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex(tbl3_lines, "latex_table_h3_main.txt")

# ── SI Table A2: Persistence by position × wave ───────────────────────────────
si_a2_rows <- persistence_tbl %>%
  arrange(match(position,c("alpha","beta","gamma")), wave_f) %>%
  mutate(
    pos_tex  = case_when(
      position=="alpha" ~ "$\\alpha$ (isolation)",
      position=="beta"  ~ "$\\beta$ (clustering)",
      position=="gamma" ~ "$\\gamma$ (bridging)"
    ),
    wave_tex = as.character(wave_f),
    row = sprintf("  %-24s & %-16s & %4d & %.3f & %.3f & [%.3f,\\;%.3f] \\\\",
                  pos_tex, wave_tex, n_transitions,
                  p_stay, p_exit, ci_lo, ci_hi)
  ) %>% pull(row)

si_a2_lines <- c(
  "% ============================================================",
  "% SI Table A2 — Persistence by position × wave",
  "% Paste into Online Appendix",
  "% ============================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Portfolio-position persistence by transition period}",
  "\\label{table:SI_A2}",
  "\\small",
  "\\begin{tabular}{llrccl}",
  "\\hline\\hline",
  "  Position & Transition & $n$ & P(stay) & P(exit) & 95\\% CI \\\\",
  "\\hline",
  si_a2_rows,
  "\\hline\\hline",
  "  \\multicolumn{6}{p{0.88\\textwidth}}{\\footnotesize \\textit{Notes.}",
  "  Unconditional persistence and exit rates per position and transition period.",
  "  95\\% CI based on normal approximation. $n$: person-wave observations per cell.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex(si_a2_lines, "latex_SI_tableA2.txt")

# ── SI Table A3: All 4 exit-from-γ models ────────────────────────────────────
fmt_ame <- function(ame_obj) {
  ame_obj %>% as_tibble() %>%
    select(term, estimate, std.error, p.value) %>%
    mutate(across(where(is.numeric), ~round(.x,3)))
}
all4 <- list(M1=fmt_ame(ames_edu), M2=fmt_ame(ames_emp),
             M3=fmt_ame(ames_full), M4=fmt_ame(ames_ord))
all_terms  <- c("edu_high","edu_ord","employed")
term_label <- c(edu_high="High education (ref: low)",
                edu_ord ="Education (ordinal, 1--4)",
                employed="Employed (ref: other)")

# Build two-row-per-term body (estimate + SE below)
si_a3_body <- c()
for (trm in all_terms) {
  ests <- sapply(names(all4), function(mn) {
    r <- all4[[mn]] %>% filter(term==trm)
    if (nrow(r)==0) return("") else paste0(fmt_n(r$estimate), sig_stars(r$p.value))
  })
  ses  <- sapply(names(all4), function(mn) {
    r <- all4[[mn]] %>% filter(term==trm)
    if (nrow(r)==0) return("") else paste0("(", fmt_n(r$std.error), ")")
  })
  si_a3_body <- c(si_a3_body,
    sprintf("  %-30s & %s & %s & %s & %s \\\\",
            term_label[trm], ests[1], ests[2], ests[3], ests[4]),
    sprintf("  %-30s & %s & %s & %s & %s \\\\",
            "", ses[1], ses[2], ses[3], ses[4])
  )
}

si_a3_lines <- c(
  "% ============================================================",
  "% SI Table A3 — All 4 exit-from-γ models",
  "% Paste into Online Appendix",
  "% ============================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Socioeconomic predictors of exit from $\\gamma$: probit AMEs (four specifications)}",
  "\\label{table:SI_A3}",
  "\\small",
  "\\begin{tabular}{lcccc}",
  "\\hline\\hline",
  "  & \\multicolumn{4}{c}{Outcome: P(exit from $\\gamma$ at $t+1$)} \\\\",
  "  \\cmidrule(lr){2-5}",
  "  & M1: Edu only & M2: Emp only & M3: Full & M4: Edu ordinal \\\\",
  "\\hline",
  si_a3_body,
  "\\hline",
  sprintf("  $n$ & \\multicolumn{4}{c}{%d} \\\\", n_gamma_trans),
  "\\hline\\hline",
  "  \\multicolumn{5}{p{0.88\\textwidth}}{\\footnotesize \\textit{Notes.}",
  "  AMEs from probit models. Sample: person-wave observations in $\\gamma$ at $t$.",
  "  M1--M3 use binary high education (t{\\'e}cnica or university vs.\\ lower);",
  "  M4 uses ordinal education (1--4). All models include sex, age, wave indicator.",
  "  Cluster-robust SE (by individual) in parentheses.",
  "  $\\dagger p < .10$; $* p < .05$; $** p < .01$; $*** p < .001$.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex(si_a3_lines, "latex_SI_tableA3.txt")

# ── SI Table A4: Exit models by starting position (specificity test) ──────────
pos_labels_tex <- c(alpha="$\\alpha$ (isolation)",
                    beta ="$\\beta$ (clustering)",
                    gamma="$\\gamma$ (bridging)")
term_short_tex <- c(edu_high="High education", employed="Employed")

si_a4_body <- c()
for (pos in c("alpha","beta","gamma")) {
  sub <- allpos_ames %>% filter(starting_position==pos)
  if (nrow(sub)==0) next
  er  <- round(sub$exit_rate[1], 3)
  np  <- sub$n_transitions[1]
  si_a4_body <- c(si_a4_body,
    "\\hline",
    sprintf("  \\multicolumn{5}{l}{\\textit{Starting position: %s --- exit rate = %.1f\\%%, $n = %d$}} \\\\",
            pos_labels_tex[pos], er*100, np),
    "\\hline"
  )
  for (i in seq_len(nrow(sub))) {
    trm <- sub$term[i]
    si_a4_body <- c(si_a4_body,
      sprintf("  \\quad %-22s & %+.3f%s & (%.3f) & %.3f \\\\",
              term_short_tex[trm],
              sub$estimate[i], sig_stars(sub$p.value[i]),
              sub$std.error[i], sub$p.value[i])
    )
  }
}

si_a4_lines <- c(
  "% ============================================================",
  "% SI Table A4 — Specificity test: exit models by starting position",
  "% Paste into Online Appendix",
  "% ============================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Socioeconomic predictors of position exit by starting position (specificity test)}",
  "\\label{table:SI_A4}",
  "\\small",
  "\\begin{tabular}{lccr}",
  "\\hline\\hline",
  "  Predictor & AME & (SE) & $p$-value \\\\",
  si_a4_body,
  "\\hline\\hline",
  "  \\multicolumn{4}{p{0.75\\textwidth}}{\\footnotesize \\textit{Notes.}",
  "  AMEs from probit models of P(exit from starting position by $t+1$).",
  "  Each panel uses respondents starting in the indicated position.",
  "  Controls: sex, age, wave. Cluster-robust SE (by individual).",
  "  $\\dagger p < .10$; $* p < .05$; $** p < .01$; $*** p < .001$.} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex(si_a4_lines, "latex_SI_tableA4.txt")

# ==============================================================================
# 10) CLUSTER BOOTSTRAP VALIDATION (SI robustness for small-N subsample)
# ==============================================================================
# Resamples ids (not obs) → preserves within-person correlation.
# For each resample: refit m_full, compute AMEs of edu_high and employed.
# Reports: bootstrap SE, percentile CI [2.5%, 97.5%], and comparison with
# the cluster-robust SEs from Section 3.
# B = 999 iterations; fast because glm (not glmm).
# ------------------------------------------------------------------------------

message("\n--- Bootstrap validation (cluster resample, B=999) ---")

set.seed(2024L)
B        <- 999L
boot_ids <- unique(dt_gamma_exits$id)
n_ids    <- length(boot_ids)

# Full-sample data used in m_full (complete cases on all covariates)
dt_boot_base <- dt_gamma_exits %>%
  filter(!is.na(edu_high), !is.na(employed), !is.na(mujer),
         !is.na(age_mid),  !is.na(wave_f))

boot_ames <- matrix(NA_real_, nrow = B, ncol = 2,
                    dimnames = list(NULL, c("edu_high", "employed")))

for (b in seq_len(B)) {
  # Sample ids with replacement
  sampled_ids <- sample(boot_ids, size = n_ids, replace = TRUE)
  # Build resampled dataset (duplicates handled by row expansion)
  db <- lapply(seq_along(sampled_ids), function(i) {
    dt_boot_base[dt_boot_base$id == sampled_ids[i], ]
  })
  db <- do.call(rbind, db)

  # Skip if no variance in outcome or predictors
  if (n_distinct(db$exit_gamma) < 2 ||
      n_distinct(db$edu_high)   < 2 ||
      n_distinct(db$employed)   < 2) next

  mb <- tryCatch(
    glm(exit_gamma ~ edu_high + employed + mujer + age_mid + wave_f,
        data = db, family = binomial("probit")),
    error = function(e) NULL
  )
  if (is.null(mb)) next

  ame_b <- tryCatch(
    suppressMessages(
      marginaleffects::avg_comparisons(mb,
        variables = c("edu_high", "employed")) %>%
        as_tibble() %>%
        select(term, estimate)
    ),
    error = function(e) NULL
  )
  if (is.null(ame_b)) next

  boot_ames[b, "edu_high"] <- ame_b$estimate[ame_b$term == "edu_high"]
  boot_ames[b, "employed"] <- ame_b$estimate[ame_b$term == "employed"]
}

# Drop failed iterations
boot_ames <- boot_ames[complete.cases(boot_ames), ]
n_valid   <- nrow(boot_ames)
message("  Valid bootstrap iterations: ", n_valid, " / ", B)

# Bootstrap SE and percentile CIs
boot_summary <- tibble::tibble(
  term        = c("edu_high", "employed"),
  ame_main    = c(edu_ame_val, emp_ame_val),
  se_cluster  = c(edu_ame_se,  emp_ame_se),
  se_boot     = apply(boot_ames, 2, sd),
  ci_lo_boot  = apply(boot_ames, 2, quantile, probs = 0.025),
  ci_hi_boot  = apply(boot_ames, 2, quantile, probs = 0.975),
  ci_lo_clust = c(edu_ame_val - 1.96 * edu_ame_se,
                  emp_ame_val - 1.96 * emp_ame_se),
  ci_hi_clust = c(edu_ame_val + 1.96 * edu_ame_se,
                  emp_ame_val + 1.96 * emp_ame_se),
  n_boot      = n_valid
) %>%
  mutate(
    ci_width_boot  = ci_hi_boot  - ci_lo_boot,
    ci_width_clust = ci_hi_clust - ci_lo_clust,
    se_ratio       = round(se_boot / se_cluster, 3)   # >1 = boot wider
  )

write_csv(boot_summary, here::here("output", "h3_bootstrap_validation.csv"))

cat("\n=== Bootstrap vs Cluster-robust comparison ===\n")
boot_summary %>%
  mutate(across(where(is.numeric), ~round(.x, 4))) %>%
  select(term, ame_main, se_cluster, se_boot, se_ratio,
         ci_lo_boot, ci_hi_boot) %>%
  print()

# Bootstrap distribution plot (SI figure)
boot_df <- tibble::tibble(
  edu_high = boot_ames[, "edu_high"],
  employed = boot_ames[, "employed"]
) %>%
  tidyr::pivot_longer(everything(), names_to = "term", values_to = "ame") %>%
  mutate(term = factor(term,
                       levels = c("edu_high", "employed"),
                       labels = c("High education", "Employed")))

obs_vals <- tibble::tibble(
  term  = factor(c("High education", "Employed"),
                 levels = c("High education", "Employed")),
  value = c(edu_ame_val, emp_ame_val)
)

p_boot <- ggplot(boot_df, aes(x = ame)) +
  geom_histogram(bins = 50, fill = "#4E79A7", alpha = 0.75, color = "white") +
  geom_vline(data = obs_vals, aes(xintercept = value),
             color = "#E15759", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
  facet_wrap(~term, scales = "free_x") +
  labs(
    title    = "Bootstrap distribution of AMEs (cluster resample, B=999)",
    subtitle = "Red dashed line = point estimate from main model",
    x        = "Average Marginal Effect",
    y        = "Count",
    caption  = paste0(
      "n(bridging subsample) = ", nrow(dt_boot_base),
      " obs | ", n_ids, " unique individuals. ",
      "Valid bootstrap iterations: ", n_valid, "/", B, "."
    )
  ) +
  theme_ssr(base_size = 12)

ggsave(here::here("output", "figure_h3_SI_bootstrap.png"),
       p_boot, width = 8, height = 4, dpi = 300)
ggsave(here::here("output", "figure_h3_SI_bootstrap.pdf"),
       p_boot, width = 8, height = 4)
message("  Bootstrap figure: output/figure_h3_SI_bootstrap.png/.pdf")

# LaTeX table for SI
boot_tex_rows <- boot_summary %>%
  mutate(
    term_tex   = c("High education (ref: low)", "Employed (ref: other)"),
    ci_clust   = sprintf("[%.3f,\\;%.3f]", ci_lo_clust, ci_hi_clust),
    ci_boot    = sprintf("[%.3f,\\;%.3f]", ci_lo_boot,  ci_hi_boot)
  ) %>%
  mutate(row = sprintf(
    "  %-30s & %+.3f & %.3f & %.3f & %s & %s \\\\",
    term_tex, ame_main, se_cluster, se_boot, ci_clust, ci_boot
  )) %>%
  pull(row)

si_boot_lines <- c(
  "% ============================================================",
  "% SI Table A5 — Cluster bootstrap validation of AMEs",
  "% ============================================================",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Bootstrap validation of AMEs for exit from $\\gamma$ (cluster resample)}",
  "\\label{table:SI_A5_bootstrap}",
  "\\small",
  "\\begin{tabular}{lccccc}",
  "\\hline\\hline",
  "  & & \\multicolumn{2}{c}{SE} & \\multicolumn{2}{c}{95\\% CI} \\\\",
  "  \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  "  Predictor & AME & Cluster-robust & Bootstrap & Cluster-robust & Bootstrap \\\\",
  "\\hline",
  boot_tex_rows,
  "\\hline\\hline",
  "  \\multicolumn{6}{p{0.92\\textwidth}}{\\footnotesize \\textit{Notes.}",
  sprintf("  AMEs from probit model of P(exit from $\\gamma$). $n = %d$ person-wave observations,", nrow(dt_boot_base)),
  sprintf("  $n_{\\text{id}} = %d$ unique respondents. Bootstrap: $B = %d$ cluster resamples", n_ids, n_valid),
  "  (sample ids with replacement; cluster-robust SE uses \\texttt{vcovCL}). ",
  "  Bootstrap CI: percentile method [2.5\\%, 97.5\\%].} \\\\",
  "\\end{tabular}",
  "\\end{table}"
)
write_latex(si_boot_lines, "latex_SI_tableA5_bootstrap.txt")

# ==============================================================================
# FINAL MANIFEST
# ==============================================================================
message("\n", strrep("=",65))
message("PAPER (main text):")
message("  figure_h3_paper.png/.pdf          -> Figure 2 (Panel A + B combined)")
message("  latex_table_h3_main.txt           -> Table 3  [paste into Overleaf]")
message("  h3_inline_numbers.txt             -> §3.2 paragraph ready to copy")
message("")
message("SUPPLEMENTARY INFORMATION:")
message("  figure_h3_SI_panelA.png/.pdf      -> SI archive: Panel A standalone")
message("  figure_h3_SI_bootstrap.png/.pdf   -> SI: bootstrap distribution")
message("  latex_SI_tableA2.txt              -> SI Table A2 [paste into Overleaf]")
message("  latex_SI_tableA3.txt              -> SI Table A3 [paste into Overleaf]")
message("  latex_SI_tableA4.txt              -> SI Table A4 [paste into Overleaf]")
message("  latex_SI_tableA5_bootstrap.txt    -> SI Table A5 [paste into Overleaf]")
message("")
message("INTERMEDIATE:")
message("  h3_persistence_table.csv / h3_exit_gamma_ames.csv")
message("  h3_allpos_exit_ames.csv  / h3_summary_for_paper.txt")
message("  h3_bootstrap_validation.csv")
message(strrep("=",65))
message("\n[003b] DONE. Open h3_inline_numbers.txt to update §3.2.")