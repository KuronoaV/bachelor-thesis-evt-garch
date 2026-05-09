# 09_conditional_evt.R
# Full-sample MLE fit of the GPD to the right tail of the standardised
# residuals, with parameter table and log-log tail-survival figure
# (Section 5.5).

source("scripts/00_setup.R")

if (!requireNamespace("rugarch",  quietly = TRUE)) install.packages("rugarch")
if (!requireNamespace("extRemes", quietly = TRUE)) install.packages("extRemes")
if (!requireNamespace("MASS",     quietly = TRUE)) install.packages("MASS")
suppressPackageStartupMessages({
  library(rugarch)
  library(extRemes)
  library(MASS)
})

spx  <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))
ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

garch_spec <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "norm"
)

fit_gpd_conditional <- function(losses, market) {
  fit_g <- ugarchfit(spec = garch_spec,
                     data = -as.numeric(coredata(losses)),
                     solver = "hybrid")
  z <- -as.numeric(residuals(fit_g, standardize = TRUE))

  u   <- quantile(z, probs = cfg$u_quantile, names = FALSE)
  fit <- fevd(z, threshold = u, type = "GP", method = "MLE")
  par <- fit$results$par
  se  <- sqrt(diag(safe_solve(fit$results$hessian)))

  scale <- par[["scale"]]; shape <- par[["shape"]]
  scale_se <- se[1];       shape_se <- se[2]

  excesses <- sort(z[z > u])
  n_total  <- length(z)
  n_exc    <- length(excesses)

  emp_surv  <- (n_exc:1) / n_total
  gpd_surv  <- (n_exc / n_total) * pmax(0, 1 + shape * (excesses - u) / scale)^(-1 / shape)
  norm_surv <- 1 - pnorm(excesses)
  t_fit     <- suppressWarnings(fitdistr(z, "t"))
  t_surv    <- 1 - pt((excesses - t_fit$estimate["m"]) / t_fit$estimate["s"],
                      df = t_fit$estimate["df"])

  surv_df <- tibble(
    Log_z      = log(excesses),
    Empirical  = log(emp_surv),
    GPD        = log(gpd_surv),
    Student_t  = log(t_surv),
    Normal     = log(norm_surv)
  ) %>%
    pivot_longer(-Log_z, names_to = "Distribution", values_to = "Log_Survival") %>%
    mutate(Distribution = factor(Distribution,
                                 levels = c("Empirical", "GPD", "Student_t", "Normal")))

  emp_color <- palette_market[[market]]
  surv_plot <- ggplot(surv_df, aes(Log_z, Log_Survival,
                                   colour = Distribution, linetype = Distribution)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c(Empirical = emp_color, GPD = "black",
                                   Student_t = "black", Normal = "black")) +
    scale_linetype_manual(values = c(Empirical = "solid", GPD = "dashed",
                                     Student_t = "dotted", Normal = "dotdash")) +
    labs(title = bquote(.(market) ~ "- Log-log tail survival (residuals)"),
         x = expression(log ~ z), y = expression(log ~ P(Z > z)))

  params <- tibble(
    Market      = market,
    Threshold_u = u,
    Exceedances = n_exc,
    Total_Obs   = n_total,
    Scale_beta  = scale,
    Scale_SE    = scale_se,
    Scale_t     = scale / scale_se,
    Scale_p     = 2 * (1 - pnorm(abs(scale / scale_se))),
    Shape_xi    = shape,
    Shape_SE    = shape_se,
    Shape_t     = shape / shape_se,
    Shape_p     = 2 * (1 - pnorm(abs(shape / shape_se)))
  )

  list(params = params, plot = surv_plot)
}

spx_cond  <- fit_gpd_conditional(spx,  "S&P 500")
ftse_cond <- fit_gpd_conditional(ftse, "FTSE MIB")

save_table(bind_rows(spx_cond$params, ftse_cond$params),
           "conditional_gpd_parameters.csv")

save_figure(spx_cond$plot | ftse_cond$plot,
            "conditional_tail_survival.svg", width = 12, height = 6)
