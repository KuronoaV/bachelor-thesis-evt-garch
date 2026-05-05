# 06_garch_estimation.R
#
# Full-sample fit of the AR(1)-GJR-GARCH(1,1) volatility model to the daily
# returns of the S&P 500 and the FTSE MIB. This is the in-sample step that
# anchors the conditional EVT pipeline of Chapter 6 and produces the parameter
# table and conditional-volatility figure for Section 5.4 of the thesis. The
# model itself is the one specified in Section 3.5: an AR(1) mean equation and
# a GJR variance equation with a leverage indicator on past negative
# innovations, fitted by Quasi-Maximum Likelihood under a Normal working
# density. Innovations are not assumed to be Normal — the QMLE step is a
# filter, not a parametric distribution claim — and the Bollerslev-Wooldridge
# robust standard errors reported here are what makes the inference valid
# under that misspecification.
#
# A point worth being explicit about: the model is fitted to the return
# series r_t = -L_t, not to the losses themselves. Section 3.5 writes the
# leverage term as gamma_1 * I[epsilon_{t-1} < 0] * epsilon_{t-1}^2 with
# gamma_1 > 0 capturing the asymmetric response of volatility to negative
# market shocks. Fitting the same specification directly on losses would flip
# the meaning of the indicator (positive losses are negative returns) and
# return a negative gamma_1, breaking the textbook reading of the parameter
# the thesis cites. Returning to losses afterwards is a sign flip on the
# fitted mean and on the standardised residuals; the conditional volatility
# is invariant to the sign of the data and reads the same in either
# convention. Subsequent scripts (07, 08, 09) inherit this convention.
#
# The figure overlays the in-sample conditional standard deviation on the
# absolute losses |L_t|, which is the cleanest visual summary of what the
# filter is doing: it tracks the slow envelope of the squared losses without
# itself being squared, and the contrast between the noisy grey series and
# the smoother coloured one is what readers should see in Section 5.4.

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
