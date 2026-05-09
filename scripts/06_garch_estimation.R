# 06_garch_estimation.R
# Full-sample AR(1)-GJR-GARCH(1,1) QMLE fit (Section 5.3). The GARCH is
# fitted to returns r_t = -L_t to keep the leverage parameter positive;
# residuals are flipped back to loss form in the scripts that consume them.

source("scripts/00_setup.R")

if (!requireNamespace("rugarch", quietly = TRUE)) install.packages("rugarch")
suppressPackageStartupMessages(library(rugarch))

spx  <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))
ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

garch_spec <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "norm"
)

fit_garch_full <- function(losses, market) {
  returns <- -as.numeric(coredata(losses))
  fit <- ugarchfit(spec = garch_spec, data = returns, solver = "hybrid")

  rc <- fit@fit$robust.matcoef
  params <- tibble(
    Market    = market,
    Parameter = rownames(rc),
    Estimate  = rc[, " Estimate"],
    Robust_SE = rc[, " Std. Error"],
    t_value   = rc[, " t value"],
    p_value   = rc[, "Pr(>|t|)"]
  )

  vol_df <- tibble(
    Date    = index(losses),
    AbsLoss = abs(as.numeric(coredata(losses))),
    Sigma   = as.numeric(sigma(fit))
  )

  vol_plot <- ggplot(vol_df, aes(Date)) +
    geom_line(aes(y = AbsLoss), colour = "grey75", linewidth = 0.3) +
    geom_line(aes(y = Sigma),   colour = palette_market[[market]], linewidth = 0.5) +
    labs(title = bquote(.(market) ~ "- " ~ "|" * L[t] * "|" ~ "vs GJR-GARCH conditional volatility"),
         x = NULL, y = "Magnitude")

  list(params = params, plot = vol_plot)
}

spx_garch  <- fit_garch_full(spx,  "S&P 500")
ftse_garch <- fit_garch_full(ftse, "FTSE MIB")

save_table(bind_rows(spx_garch$params, ftse_garch$params),
           "garch_parameters.csv")

save_figure(spx_garch$plot / ftse_garch$plot,
            "garch_conditional_volatility.svg", width = 10, height = 8)
