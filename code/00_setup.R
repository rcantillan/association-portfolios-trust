# ==============================================================================
# 00_setup.R — Replication package setup
# ------------------------------------------------------------------------------
# - Sets project root via {here}
# - Loads/installs required packages
# - Defines common paths and small helpers
# ==============================================================================

rm(list = ls())

# ---- Packages ----
req_pkgs <- c(
  "here", "data.table", "tidyverse", "LMest",
  "lme4", "marginaleffects", "broom", "broom.mixed",
  "nnet", "patchwork"
)

installed <- rownames(installed.packages())
to_install <- setdiff(req_pkgs, installed)
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

invisible(lapply(req_pkgs, library, character.only = TRUE))

# ---- Paths ----
ROOT   <- here::here()
DATA   <- here::here("data")
CODE   <- here::here("code")
OUT    <- here::here("output")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

message("Project root: ", ROOT)
message("Data dir:      ", DATA)
message("Output dir:    ", OUT)

# ---- Helpers ----
# Safer ggsave wrapper
ggsave_safe <- function(filename, plot, width = 8, height = 5, dpi = 300) {
  ggplot2::ggsave(filename = here::here("output", filename),
                  plot = plot, width = width, height = height, dpi = dpi)
}

# Convert to 0/1 robustly
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x > 0 & !is.na(x)))
  x2 <- tolower(trimws(as.character(x)))
  as.integer(x2 %in% c("1","yes","y","true","t","si","sí"))
}

# Simple check
stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing required file(s):\n- ", paste(missing, collapse = "\n- "))
  }
}
