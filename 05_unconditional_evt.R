# 05_unconditional_evt.R
#
# Full-sample fit of the Generalised Pareto Distribution to the right tail of
# the raw daily losses, applied independently to the S&P 500 and the FTSE MIB.
# This is the in-sample counterpart of the Unconditional EVT model evaluated
# out-of-sample in Chapter 6 (rolling-window n = 2500). The point of running
# the fit on the entire history once is to populate the in-sample diagnostic
# section of Chapter 5 and to give the reader a reference value for the
# shape parameter against which the rolling estimates of script 16 can be
# read.
#
# The threshold is fixed at the 90th percentile of the loss distribution,
# the choice defended in Section 4.4 and supported by the diagnostics of
# script 04. Maximum-likelihood estimates of the scale and shape parameters
# are reported with robust Wald standard errors derived from the inverse
# Hessian, together with t-statistics and p-values.
#
# The accompanying graphical diagnostic is the GPD QQ plot, drawn on the
# excesses Y_j = L_j - u for L_j > u. Departure of the empirical excesses
# from the 45-degree line in the upper right corner — the regime where the
# GPD makes its strongest predictions — is the diagnostic that matters most
# for risk forecasting. The numeric outputs are saved as plain numbers (not
# pre-formatted strings) so that Chapter 5 can format them consistently with
# the rest of the table set.

source("scripts/00_setup.R")

if (!requireNamespace("extRemes", quietly = TRUE)) install.packages("extRemes")
suppressPackageStartupMessages(library(extRemes))

spx  <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))
ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

fit_gpd_full <- function(x, market) {
  x <- as.numeric(coredata(x))
  u <- quantile(x, probs = cfg$u_quantile, names = FALSE)

  fit <- fevd(x, threshold = u, type = "GP", method = "MLE")
  par <- fit$results$par
  se  <- sqrt(diag(safe_solve(fit$results$hessian)))

  scale <- par[["scale"]]; shape <- par[["shape"]]
  scale_se <- se[1];       shape_se <- se[2]

  excesses <- sort(x[x > u] - u)
  n_exc    <- length(excesses)
  p_seq    <- (seq_len(n_exc) - 0.5) / n_exc
  theo_q   <- (scale / shape) * ((1 - p_seq)^(-shape) - 1)

  qq <- ggplot(tibble(Theoretical = theo_q, Empirical = excesses),
               aes(Theoretical, Empirical)) +
    geom_point(colour = palette_market[[market]], alpha = 0.5, shape = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.7) +
    labs(title = bquote(.(market) ~ "- GPD QQ (raw losses)"),
         x = "Theoretical GPD quantiles", y = "Empirical excesses")

  params <- tibble(
    Market      = market,
    Threshold_u = u,
    Exceedances = n_exc,
    Total_Obs   = length(x),
    Scale_beta  = scale,
    Scale_SE    = scale_se,
    Scale_t     = scale / scale_se,
    Scale_p     = 2 * (1 - pnorm(abs(scale / scale_se))),
    Shape_xi    = shape,
    Shape_SE    = shape_se,
    Shape_t     = shape / shape_se,
    Shape_p     = 2 * (1 - pnorm(abs(shape / shape_se)))
  )

  list(params = params, plot = qq)
}

spx_fit  <- fit_gpd_full(spx,  "S&P 500")
ftse_fit <- fit_gpd_full(ftse, "FTSE MIB")

save_table(bind_rows(spx_fit$params, ftse_fit$params),
           "unconditional_gpd_parameters.csv")

save_figure(spx_fit$plot | ftse_fit$plot,
            "unconditional_gpd_qq.svg", width = 10, height = 5)
