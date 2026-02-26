# ==============================================================================
# 03_H4_connectivity_figures_themeSSR.R — H4 figures with your theme_ssr()
# ------------------------------------------------------------------------------
# Outputs:
#  - H4_prev_gamma_soft.(png|pdf)
#  - H4_mean_jacc.(png|pdf)
#  - H4_panel_gamma_jacc.(png|pdf)  [recommended for SI; can also be main-text]
# ==============================================================================

source(here::here("code","00_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(readr)
  library(patchwork)  # for combined panels
})

dt <- readRDS(here::here("data","dt_states.rds"))
dt <- as.data.table(dt)

domains <- c("nhg","religious","political","union","professional","charity","sport","student")

# Guardrails
miss <- setdiff(c("id","ola", domains, "p_gamma", "position"), names(dt))
if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse=", "))

# ---- Compute H4 metrics (NA-robust) ----
h4 <- dt[, {
  X <- as.matrix(.SD[, ..domains])
  X[is.na(X)] <- 0
  storage.mode(X) <- "numeric"

  P <- t(X) %*% X
  diagv <- diag(P)

  # mean pairwise Jaccard overlap across domain pairs
  J <- matrix(NA_real_, nrow = length(domains), ncol = length(domains))
  for (a in seq_along(domains)) for (b in seq_along(domains)) {
    if (a == b) next
    union_ab <- diagv[a] + diagv[b] - P[a, b]
    J[a, b] <- ifelse(union_ab > 0, P[a, b] / union_ab, NA_real_)
  }
  mean_jacc <- mean(J[upper.tri(J)], na.rm = TRUE)

  # simple multi-domain prevalence
  dom_count <- rowSums(X)
  p_two_plus <- mean(dom_count >= 2)

  prev_gamma_soft <- mean(p_gamma, na.rm = TRUE)
  prev_gamma_hard <- mean(position == "gamma", na.rm = TRUE)

  .(mean_jacc = mean_jacc,
    p_two_plus = p_two_plus,
    prev_gamma_soft = prev_gamma_soft,
    prev_gamma_hard = prev_gamma_hard)
}, by = ola, .SDcols = domains][order(ola)]

dir.create(here::here("output"), showWarnings = FALSE, recursive = TRUE)
write_csv(h4, here::here("output","H4_connectivity_by_wave_jaccard.csv"))

# ---- Figure 1: Bridging prevalence (soft) ----
p_gamma <- ggplot(h4, aes(x = ola, y = prev_gamma_soft, group = 1)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(h4$prev_gamma_soft, na.rm = TRUE) * 1.15)
  ) +
  labs(
    x = "Wave",
    y = "Pr(bridging) [mean posterior p_gamma]",
    title = "Bridging prevalence by wave (soft)"
  ) +
  theme_ssr()

# ---- Figure 2: Connectivity (mean Jaccard) ----
p_jacc <- ggplot(h4, aes(x = ola, y = mean_jacc, group = 1)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(
    x = "Wave",
    y = "Mean pairwise Jaccard overlap (domains)",
    title = "Portfolio-implied cross-domain connectivity by wave"
  ) +
  theme_ssr()

# ---- Combined panel (recommended for SI; can also be main) ----
p_panel <- (p_gamma / p_jacc) +
  plot_annotation(
    title = "Portfolio dynamics and portfolio-implied cross-domain connectivity",
    theme = theme_ssr()
  )

# ---- Save (your helpers + explicit PDFs) ----
ggsave_safe("H4_prev_gamma_soft.png", p_gamma, width = 6.8, height = 4.4)
ggsave_safe("H4_mean_jacc.png",       p_jacc,  width = 6.8, height = 4.4)
ggsave_safe("H4_panel_gamma_jacc.png",p_panel, width = 7.2, height = 8.6)

ggsave(here::here("output","H4_prev_gamma_soft.pdf"), p_gamma, width = 6.8, height = 4.4)
ggsave(here::here("output","H4_mean_jacc.pdf"),       p_jacc,  width = 6.8, height = 4.4)
ggsave(here::here("output","H4_panel_gamma_jacc.pdf"),p_panel, width = 7.2, height = 8.6)

message("DONE: H4 figures saved with theme_ssr() in output/")
message("CSV: output/H4_connectivity_by_wave_jaccard.csv")