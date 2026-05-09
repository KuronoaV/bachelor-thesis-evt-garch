# 07_residual_diagnostics.R
# Ljung-Box on z_t and z_t^2, Engle's ARCH-LM at lags 10/15/20, and ACF
# plots of the standardised GARCH residuals (Section 5.4).

source("scripts/00_setup.R")

if (!requireNamespace("rugarch", quietly = TRUE)) install.packages("rugarch")
if (!requireNamespace("FinTS",   quietly = TRUE)) install.packages("FinTS")
suppressPackageStartupMessages({
  library(rugarch)
  library(FinTS)
})

spx  <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))
ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

garch_spec <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "norm"
)

residual_tests <- function(z, market, lags = c(10, 15, 20)) {
  bind_rows(lapply(lags, function(l) {
    lb_z  <- Box.test(z,   lag = l, type = "Ljung-Box")
    lb_z2 <- Box.test(z^2, lag = l, type = "Ljung-Box")
    arch  <- ArchTest(z, lags = l)
    tibble(
      Market       = market,
      Lag          = l,
      LB_Raw_Stat  = unname(lb_z$statistic),
      LB_Raw_p     = lb_z$p.value,
      LB_Sq_Stat   = unname(lb_z2$statistic),
      LB_Sq_p      = lb_z2$p.value,
      ARCH_LM_Stat = unname(arch$statistic),
      ARCH_LM_p    = arch$p.value
    )
  }))
}

acf_frame <- function(z, lag_max = 40) {
  out <- acf(z, plot = FALSE, lag.max = lag_max)
  tibble(lag = as.numeric(out$lag[-1]), acf = as.numeric(out$acf[-1]))
}

acf_plot <- function(df, n_obs, title, market) {
  ci <- 1.96 / sqrt(n_obs)
  ggplot(df, aes(lag, acf)) +
    geom_segment(aes(xend = lag, yend = 0),
                 linewidth = 0.3, colour = palette_market[[market]]) +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    coord_cartesian(ylim = c(-0.1, 0.1)) +
    labs(title = title, x = "Lag", y = "Autocorrelation")
}

run_diagnostics <- function(losses, market) {
  fit <- ugarchfit(spec = garch_spec,
                   data = -as.numeric(coredata(losses)),
                   solver = "hybrid")
  z <- as.numeric(residuals(fit, standardize = TRUE))

  list(
    table = residual_tests(z, market),
    p_z   = acf_plot(acf_frame(z),    length(z), bquote(.(market) ~ "- ACF of" ~ z[t]),    market),
    p_z2  = acf_plot(acf_frame(z^2),  length(z), bquote(.(market) ~ "- ACF of" ~ z[t]^2), market)
  )
}

spx_diag  <- run_diagnostics(spx,  "S&P 500")
ftse_diag <- run_diagnostics(ftse, "FTSE MIB")

save_table(bind_rows(spx_diag$table, ftse_diag$table),
           "residual_diagnostics.csv")

save_figure(
  (spx_diag$p_z  | ftse_diag$p_z) /
  (spx_diag$p_z2 | ftse_diag$p_z2),
  "residual_acf_plots.svg", width = 12, height = 8
)
