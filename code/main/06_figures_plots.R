# ── 06_figures_plots.R — Main manuscript figures ──────────────────────────────
# Requires in memory:
#   prob_df: data.frame with columns item + State_1..State_3 (P(Y=1|state))
#   Pi_avg:  3x3 matrix of average transition probabilities
#   K_BASELINE (3), MOD_TRANS (1), ggsave_safe()
#
# Outputs:
#   output/fig_profiles_clean.(png/pdf)
#   output/fig_transition_ggraph.(png/pdf)
#   output/fig_transition_heatmap.(png/pdf)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(scales)
  library(tidygraph)
  library(ggraph)
  library(grid)
  library(here)
})

# ── Theme (from 00_setup.R via source) ──────────────────────────────────────
# theme_ssr() is already loaded; this script just calls it.

# ── Labels ───────────────────────────────────────────────────────────────────
domain_order <- c("nhg","religious","political","union","professional","charity","sport","student")

domain_labels2 <- c(
  nhg          = "Neighborhood orgs.",
  religious    = "Religious orgs.",
  political    = "Political parties",
  union        = "Labor unions",
  professional = "Professional assoc.",
  charity      = "Charitable orgs.",
  sport        = "Sports clubs",
  student      = "Student orgs."
)

# State colors: use STATE_PALETTE_LEGACY from 00_setup.R (NPG palette)
state_key <- tibble(
  state_raw   = c("State_1", "State_2", "State_3"),
  state_long  = c("α (isolation)", "β (clustering)", "γ (bridging)"),
  state_short = c("α", "β", "γ"),
  # NPG colors: α=navy, β=teal, γ=red
  col = c("#3C5488", "#4DBBD5", "#E64B35")
)

state_palette <- setNames(state_key$col, state_key$state_long)
state_levels  <- state_key$state_long
states_short  <- state_key$state_short

# Facet label: "α  Isolation  [55%]" style — added as factor label
STATE_SHARE <- tryCatch({
  # initial distribution from Pi_avg (row averages ≈ steady-state shares)
  tibble(
    state_long  = c("α (isolation)", "β (clustering)", "γ (bridging)"),
    share_label = c(
      sprintf("α — Isolation  [%d%%]",  round(mean(Pi_avg[, "alpha"  ]) * 100)),
      sprintf("β — Clustering [%d%%]",  round(mean(Pi_avg[, "beta"   ]) * 100)),
      sprintf("γ — Bridging   [%d%%]",  round(mean(Pi_avg[, "gamma"  ]) * 100))
    )
  )
}, error = function(e) {
  tibble(
    state_long  = c("α (isolation)", "β (clustering)", "γ (bridging)"),
    share_label = c("α — Isolation", "β — Clustering", "γ — Bridging")
  )
})

# ── Figure 1 — Membership profiles (lollipop) ────────────────────────────────
prob_long <- prob_df |>
  pivot_longer(starts_with("State_"), names_to = "state_raw", values_to = "prob") |>
  left_join(state_key, by = "state_raw") |>
  left_join(STATE_SHARE, by = "state_long") |>
  mutate(
    state_long   = factor(state_long, levels = state_levels),
    share_label  = factor(share_label, levels = STATE_SHARE$share_label),
    item         = factor(item, levels = rev(domain_order)),
    item_lab     = recode(as.character(item), !!!domain_labels2),
    pct_label    = scales::percent(prob, accuracy = 1)
  )

p_profiles_clean <- ggplot(
  prob_long,
  aes(x = item_lab, y = prob, color = state_long)
) +
  # stem from 0 to point
  geom_segment(
    aes(xend = item_lab, yend = 0),
    color    = "grey82",
    linewidth = 0.55
  ) +
  # dot
  geom_point(size = 3.2, show.legend = FALSE) +
  # percentage label to the right of dot
  geom_text(
    aes(label = pct_label),
    hjust  = -0.32,
    size   = 3.0,
    color  = "grey35",
    family = "sans"
  ) +
  coord_flip() +
  facet_wrap(~share_label, ncol = 1) +
  scale_color_manual(values = state_palette) +
  scale_y_continuous(
    limits = c(0, 1.12),
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Pr(member | state)") +
  theme_ssr(base_size = 13) +
  theme(
    legend.position    = "none",
    strip.text         = element_text(size = 12.5, color = "grey15"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey93", linewidth = 0.3),
    axis.text.y        = element_text(size = 11.5, color = "grey25"),
    panel.spacing      = unit(1.1, "lines")
  )

ggsave_safe("fig_profiles_clean.png", p_profiles_clean, width = 7.5, height = 9.5)
ggsave(here::here("output", "fig_profiles_clean.pdf"), p_profiles_clean, width = 7.5, height = 9.5)

# ── Figure 2 — Transition network ─────────────────────────────────────────────
# 9 directed transitions from Pi_avg; arrows trimmed to node ring (NODE_R).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggforce)
  library(grid)
  library(here)
})

# Tunable knobs
CURV_OUT  <- 1.25
CURV_IN   <- 0.65

LAB_T     <- 0.55
LAB_NOFF  <- 0.06

# Arrowheads as curved mini-beziers near the end
AR_T0     <- 0.90
AR_T1     <- 0.995
ARROW_SZ  <- unit(6, "mm")

# Node ring radius (THIS is what makes arrow tips touch the ring)
LOOP_R    <- 0.30
NODE_R    <- LOOP_R       # start/end trimming radius; set = LOOP_R by default

LOOP_ANG  <- c("α" = 110, "β" = 160, "γ" =  20)

# Width variation
W_MIN <- 0.4
W_MAX <- 7.0
W_POW <- 0.70

A_MIN <- 0.20
A_MAX <- 0.95

LAB_DIGITS <- 3
fmt <- function(x) sprintf(paste0("%.", LAB_DIGITS, "f"), x)

# State labels and colors
states_short <- c("α","β","γ")
state_levels <- c("α (isolation)","β (clustering)","γ (bridging)")
# NPG colors (match STATE_PALETTE_LEGACY from 00_setup.R)
state_palette <- c(
  "α (isolation)"  = "#3C5488",   # dark navy
  "β (clustering)" = "#4DBBD5",   # bright teal
  "γ (bridging)"   = "#E64B35"    # warm red
)

# Nodes in triangle layout
nodes <- tibble(
  name = states_short,
  state_long = state_levels,
  x = c(0, -2,  2),
  y = c(2,  0,  0)
)

# 1) Extract all 9 transitions from Pi_avg
edges <- expand.grid(from_i = 1:3, to_i = 1:3) %>%
  as_tibble() %>%
  mutate(
    p = as.numeric(Pi_avg[cbind(from_i, to_i)]),
    p = if_else(is.na(p), 0, p),
    from = states_short[from_i],
    to   = states_short[to_i],
    is_loop = (from_i == to_i),
    label = fmt(p),
    lw = W_MIN + (W_MAX - W_MIN) * (p ^ W_POW),
    a  = pmin(A_MAX, pmax(A_MIN, p))
  ) %>%
  select(from_i, to_i, from, to, p, label, is_loop, lw, a)

stopifnot(nrow(edges) == 9)

# 2) Join node coordinates and compute geometry
E <- edges %>%
  left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>% rename(x1 = x, y1 = y) %>%
  left_join(nodes %>% select(name, x, y), by = c("to"   = "name")) %>% rename(x2 = x, y2 = y) %>%
  mutate(
    dx = x2 - x1, dy = y2 - y1,
    denom = pmax(1e-9, sqrt(dx^2 + dy^2)),
    nx = -dy / denom,
    ny =  dx / denom,

    # unordered pair id (numeric)
    pair = paste0(pmin(from_i, to_i), pmax(from_i, to_i)),
    # direction within pair (+1 low→high, -1 high→low, 0 for loops)
    dir  = if_else(from_i < to_i, 1, if_else(from_i > to_i, -1, 0)),

    # which side is "outside" for each unordered pair (tuned to your triangle)
    pair_side = case_when(
      pair == "12" ~ -1,  # α-β : outside to the left
      pair == "13" ~  1,  # α-γ : outside up/right
      pair == "23" ~ -1,  # β-γ : outside down
      TRUE ~ 1
    ),

    # separates directions
    curv_sign = dir * pair_side,
    curv_mag  = if_else(curv_sign == 1, CURV_OUT, CURV_IN)
  )

# 3) Non-loop edges: trim endpoints to touch the node ring
E_nl <- E %>%
  filter(!is_loop) %>%
  mutate(
    # chord unit direction
    ux = dx / denom,
    uy = dy / denom,

    # TRIM: move start/end from centers to the ring boundary
    xs = x1 + NODE_R * ux,
    ys = y1 + NODE_R * uy,
    xe = x2 - NODE_R * ux,
    ye = y2 - NODE_R * uy,

    # midpoint of trimmed segment
    mx2 = (xs + xe) / 2,
    my2 = (ys + ye) / 2,

    # bezier control point (use same normal, but trimmed endpoints)
    cx = mx2 + curv_sign * curv_mag * nx,
    cy = my2 + curv_sign * curv_mag * ny,

    gid = row_number()
  )

# Quadratic bezier point at t (using trimmed endpoints)
bez_xy <- function(xs,ys,cx,cy,xe,ye,t){
  tibble(
    x = (1-t)^2 * xs + 2*(1-t)*t*cx + t^2*xe,
    y = (1-t)^2 * ys + 2*(1-t)*t*cy + t^2*ye
  )
}

# Full bezier paths (3 points)
bezier_df <- bind_rows(
  E_nl %>% transmute(gid, x = xs, y = ys, lw, a),
  E_nl %>% transmute(gid, x = cx, y = cy, lw, a),
  E_nl %>% transmute(gid, x = xe, y = ye, lw, a)
)

# Arrowhead curved mini-bezier near the end (t0, tm, t1)
t0 <- AR_T0
t1 <- AR_T1
tm <- (t0 + t1) / 2

arrow_bezier_df <- bind_rows(
  E_nl %>% rowwise() %>% mutate(p0 = list(bez_xy(xs,ys,cx,cy,xe,ye,t0))) %>% unnest_wider(p0) %>%
    transmute(gid, x = x, y = y, lw = lw, step = 1),
  E_nl %>% rowwise() %>% mutate(pm = list(bez_xy(xs,ys,cx,cy,xe,ye,tm))) %>% unnest_wider(pm) %>%
    transmute(gid, x = x, y = y, lw = lw, step = 2),
  E_nl %>% rowwise() %>% mutate(p1 = list(bez_xy(xs,ys,cx,cy,xe,ye,t1))) %>% unnest_wider(p1) %>%
    transmute(gid, x = x, y = y, lw = lw, step = 3)
) %>%
  arrange(gid, step)

# Labels on directed non-loop edges (embedded)
label_nl <- E_nl %>%
  rowwise() %>%
  mutate(
    pL = list(bez_xy(xs,ys,cx,cy,xe,ye,LAB_T)),
    lx = pL$x + curv_sign * LAB_NOFF * nx,
    ly = pL$y + curv_sign * LAB_NOFF * ny
  ) %>%
  ungroup() %>%
  select(lx, ly, label)

# 4) Loops: circles + labels
E_lp <- E %>% filter(is_loop)

loop_df <- E_lp %>%
  transmute(
    x0 = x1, y0 = y1, r = LOOP_R,
    lw, a, label,
    node = from
  )

loop_labels <- loop_df %>%
  rowwise() %>%
  mutate(
    ang = (pi/180) * LOOP_ANG[node],
    lx = x0 + r * cos(ang),
    ly = y0 + r * sin(ang)
  ) %>%
  ungroup() %>%
  select(lx, ly, label)

# 5) Build plot
p_trans <- ggplot() +

  # Curved edges (base)
  ggforce::geom_bezier(
    data = bezier_df,
    aes(x = x, y = y, group = gid, linewidth = lw, alpha = a),
    colour = "grey60",
    lineend = "round",
    show.legend = FALSE
  ) +

  # Curved arrowhead segment (draw on top, opaque for crisp heads)
  ggforce::geom_bezier(
    data = arrow_bezier_df,
    aes(x = x, y = y, group = gid, linewidth = lw),
    colour = "grey60",
    alpha  = 2,
    arrow  = grid::arrow(length = ARROW_SZ, type = "closed"),
    lineend = "round",
    show.legend = FALSE
  ) +

  # Loop circles
  ggforce::geom_circle(
    data = loop_df,
    aes(x0 = x0, y0 = y0, r = r, linewidth = lw, alpha = a),
    colour = "grey60",
    show.legend = FALSE
  ) +

  # Edge labels — no border, white fill blends with background
  geom_label(
    data = label_nl,
    aes(x = lx, y = ly, label = label),
    fill        = "white",
    color       = "grey25",
    label.size  = 0,          # no border
    label.padding = unit(0.08, "lines"),
    size        = 3.5,
    family      = "sans",
    show.legend = FALSE
  ) +
  geom_label(
    data = loop_labels,
    aes(x = lx, y = ly, label = label),
    fill        = "white",
    color       = "grey25",
    label.size  = 0,
    label.padding = unit(0.08, "lines"),
    size        = 3.5,
    family      = "sans",
    show.legend = FALSE
  ) +

  # Nodes — clean circles, NPG fill, thin stroke
  geom_point(
    data = nodes,
    aes(x = x, y = y, fill = state_long),
    shape = 21, size = 20, stroke = 0.5, colour = "white",
    show.legend = FALSE
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = name),
    size = 7, fontface = "bold", color = "white"
  ) +

  scale_fill_manual(values = state_palette) +
  scale_linewidth_identity() +
  scale_alpha_identity() +
  coord_equal(clip = "off") +
  theme_void(base_size = 14) +
  theme(plot.margin = margin(12, 30, 12, 30)) +
  labs()

# Save
ggsave(here::here("output","fig_transition_ggraph.pdf"), p_trans, width = 7.8, height = 5.6)
ggsave(here::here("output","fig_transition_ggraph.png"), p_trans, width = 7.8, height = 5.6, dpi = 300)

p_trans

# ── Figure 3 — Transition heatmap (SI) ───────────────────────────────────────
Pi_heat <- as.data.frame(as.table(Pi_avg))
colnames(Pi_heat) <- c("from_i","to_i","p")

Pi_heat <- Pi_heat %>%
  mutate(
    from = states_short[as.integer(from_i)],
    to   = states_short[as.integer(to_i)],
    lab  = sprintf("%.2f", p)
  )

p_trans_heat <- ggplot(Pi_heat, aes(x = to, y = from, fill = p)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = lab), size = 6.2) +
  scale_fill_gradient(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  labs(
    title = paste0("Transition probabilities (average) (K=", K_BASELINE, ", mod=", MOD_TRANS, ")"),
    x = "To state",
    y = "From state"
  ) +
  theme_ssr_big() +
  theme(legend.position = "right")

ggsave_safe("fig_transition_heatmap.png", p_trans_heat, width = 8.2, height = 6.2)
ggsave(here::here("output","fig_transition_heatmap.pdf"), p_trans_heat, width = 8.2, height = 6.2)