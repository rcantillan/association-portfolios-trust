# ==============================================================================
# 00_setup.R — Global setup (packages, paths, helpers)
# ------------------------------------------------------------------------------
# All scripts assume they are run from repository root.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(data.table)
})

# ---- Paths ----
DIR_DATA <- here::here("data")
DIR_OUT  <- here::here("output")
if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)

# ---- Required file check ----
stop_if_missing <- function(paths) {
  miss <- paths[!file.exists(paths)]
  if (length(miss) > 0) {
    stop("Missing required files:\n- ", paste(miss, collapse = "\n- "))
  }
}

# ---- Safe ggsave into output/ ----
ggsave_safe <- function(filename, plot, width = 7, height = 5, dpi = 300) {
  ggplot2::ggsave(
    filename = here::here("output", filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

# ---- Missing codes common in ELSOC ----
to_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- dplyr::na_if(x, -888)
  x <- dplyr::na_if(x, -999)
  x <- dplyr::na_if(x, 777777)
  x <- dplyr::na_if(x, 999999)
  x
}

# ---- Robust 0/1 conversion ----
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.character(x)) {
    x <- trimws(tolower(x))
    return(dplyr::case_when(
      x %in% c("si","sí","yes","y","1","true") ~ 1L,
      x %in% c("no","0","false") ~ 0L,
      TRUE ~ NA_integer_
    ))
  }
  x <- to_na(x)
  dplyr::case_when(
    is.na(x) ~ NA_integer_,
    x %in% c(0, 1) ~ as.integer(x),
    TRUE ~ as.integer(x)
  )
}

# ---- Global knobs ----
# Choose neighborhood trust item here and it propagates to all scripts.
TRUST_NH_VAR <- "c03"   # <-- change if needed (e.g., "c04")

# Membership items in the zip-based ELSOC longitudinal file:
MEMBER_ITEMS <- paste0("c12_0", 1:8)

# Analysis waves used in original workflow:
WAVES_RAW <- c(1, 3, 6)  # will be recoded to 1,2,3


theme_ssr_big <- function(base_size = 16, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 3),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = base_size),
      axis.text = ggplot2::element_text(size = base_size),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size),
      legend.position = "top"
    )
}
