# 00_setup.R
# Loads packages, declares paths, the cfg list (confidence levels, window
# lengths, crisis dates, MC seed), the colour palette and the ggplot theme.
# Sourced first by every other script.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(xts)
  library(zoo)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

paths <- list(
  root      = ".",
  data_raw  = "data/raw",
  data_proc = "data/processed",
  tab       = "outputs/tables",
  fig       = "outputs/figures"
)

for (p in paths[c("data_proc", "tab", "fig")]) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

cfg <- list(
  alpha_var   = 0.99,
  alpha_es    = 0.975,
  u_quantile  = 0.90,
  windows = list(
    gauss_hs = 500,
    uevt     = 2500,
    uevt_min = 1500,
    cevt     = 1000
  ),
  sample = list(
    spx_start = as.Date("1990-01-01"),
    end       = as.Date("2026-04-21")
  ),
  crises = list(
    Dotcom = list(
      label = "Dot-com Bubble",
      slug  = "dotcom",
      start = as.Date("2000-03-24"),
      end   = as.Date("2002-10-09")
    ),
    GFC = list(
      label = "Global Financial Crisis",
      slug  = "gfc",
      start = as.Date("2007-10-09"),
      end   = as.Date("2009-03-09")
    ),
    Covid = list(
      label = "COVID-19 Crash",
      slug  = "covid",
      start = as.Date("2020-02-19"),
      end   = as.Date("2020-03-23"))
  ),
  mc = list(n_sim = 10000, seed = 20260424)
)

set.seed(cfg$mc$seed)

palette_market <- c("S&P 500" = "#1f77b4", "FTSE MIB" = "#ff7f0e")
palette_model  <- c(
  "Gaussian Parametric"   = "#4C72B0",
  "Historical Simulation" = "#55A868",
  "Unconditional EVT"     = "#C44E52",
  "Conditional EVT"       = "#8172B2"
)
model_labels <- c(
  Gauss = "Gaussian Parametric",
  HS    = "Historical Simulation",
  UEVT  = "Unconditional EVT",
  CEVT  = "Conditional EVT"
)

theme_thesis <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      text             = element_text(family = "Times New Roman"),
      plot.title       = element_text(face = "plain", hjust = 0.5),
      axis.title       = element_text(face = "italic"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text       = element_text(face = "plain", size = base_size - 2),
      legend.position  = "bottom"
    )
}
theme_set(theme_thesis())

save_table <- function(x, name) {
  write.csv(x, file.path(paths$tab, name), row.names = FALSE)
  invisible(x)
}

save_figure <- function(plot, name, width = 10, height = 6) {
  ggsave(file.path(paths$fig, name), plot = plot, width = width, height = height)
  invisible(plot)
}

# Hessian inversion can fail at the boundary of the GPD parameter space; this
# returns NA standard errors instead of crashing the script.
safe_solve <- function(H) {
  out <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(out)) matrix(NA_real_, nrow(H), ncol(H)) else out
}

crisis_subset <- function(x, crisis) {
  x[paste0(format(crisis$start), "/", format(crisis$end))]
}
