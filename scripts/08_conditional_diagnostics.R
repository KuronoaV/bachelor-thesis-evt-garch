# 08_conditional_diagnostics.R
# POT threshold diagnostics (mean-excess, shape and modified-scale stability)
# computed on the standardised GARCH residuals (Section 5.5).

source("scripts/00_setup.R")

if (!requireNamespace("rugarch",  quietly = TRUE)) install.packages("rugarch")
if (!requireNamespace("extRemes", quietly = TRUE)) install.packages("extRemes")
suppressPackageStartupMessages({
  library(rugarch)
  library(extRemes)
})

spx  <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))
ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

garch_spec <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "norm"
)

loss_innovations <- function(losses) {
  fit <- ugarchfit(spec = garch_spec,
                   data = -as.numeric(coredata(losses)),
                   solver = "hybrid")
  -as.numeric(residuals(fit, standardize = TRUE))
}

mef_grid <- function(z, qs = seq(0.50, 0.995, by = 0.005)) {
  us <- quantile(z, probs = qs, names = FALSE)
  tibble(
    Quantile  = qs,
    Threshold = us,
    MEF       = sapply(us, function(u) mean(z[z > u] - u))
  )
}

# Wald CIs straight off the Hessian. Modified-scale CI uses the delta method
# so that Cov(beta, xi) is propagated correctly:
#   Var(beta*) = Var(beta) + u^2 Var(xi) - 2u Cov(beta, xi)
stability_grid <- function(z, qs = seq(0.50, 0.985, by = 0.005), min_exc = 25) {
  us <- quantile(z, probs = qs, names = FALSE)

  rows <- lapply(seq_along(us), function(i) {
    u <- us[i]
    if (sum(z > u) < min_exc) return(NULL)
    fit <- tryCatch(
      suppressWarnings(fevd(z, threshold = u, type = "GP", method = "MLE")),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    par   <- fit$results$par
    H_inv <- safe_solve(fit$results$hessian)
    scale <- par[["scale"]]; shape <- par[["shape"]]
    shape_se     <- sqrt(H_inv[2, 2])
    beta_star    <- scale - shape * u
    beta_star_var <- H_inv[1, 1] + u^2 * H_inv[2, 2] - 2 * u * H_inv[1, 2]
    beta_star_se <- sqrt(max(0, beta_star_var))

    tibble(
      Quantile     = qs[i],
      Threshold    = u,
      Shape        = shape,
      Shape_Lo     = shape - 1.96 * shape_se,
      Shape_Hi     = shape + 1.96 * shape_se,
      Adj_Scale    = beta_star,
      Adj_Scale_Lo = beta_star - 1.96 * beta_star_se,
      Adj_Scale_Hi = beta_star + 1.96 * beta_star_se
    )
  })
  bind_rows(rows)
}

mef_plot <- function(df, u_ref, market) {
  ggplot(df, aes(Threshold, MEF)) +
    geom_point(size = 1.2, colour = palette_market[[market]]) +
    geom_vline(xintercept = u_ref, linetype = "dashed", colour = "red") +
    annotate("text", x = u_ref, y = max(df$MEF),
             label = "u == q[0.90]", parse = TRUE,
             colour = "red", hjust = -0.15, size = 3) +
    labs(title = bquote(.(market) ~ "- Mean excess function (residuals)"),
         x = "Threshold (u)", y = "Mean excess")
}

stability_plot <- function(df, u_ref, market, value, lo, hi, kind) {
  if (kind == "shape") {
    title_expr <- bquote(.(market) ~ "-" ~ xi ~ "stability (residuals)")
    y_expr     <- expression(xi)
  } else {
    title_expr <- bquote(.(market) ~ "-" ~ beta * "*" ~ "stability (residuals)")
    y_expr     <- expression(beta * "*")
  }
  ggplot(df, aes(Threshold, .data[[value]])) +
    geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]]),
                fill = palette_market[[market]], alpha = 0.15) +
    geom_line(linewidth = 0.5, colour = palette_market[[market]]) +
    geom_vline(xintercept = u_ref, linetype = "dashed", colour = "red") +
    labs(title = title_expr, x = "Threshold (u)", y = y_expr)
}

build_panel <- function(losses, market) {
  z      <- loss_innovations(losses)
  u_ref  <- quantile(z, probs = cfg$u_quantile, names = FALSE)
  mef    <- mef_grid(z)
  stab   <- stability_grid(z)

  p_mef   <- mef_plot(mef, u_ref, market)
  p_shape <- stability_plot(stab, u_ref, market, "Shape",
                            "Shape_Lo", "Shape_Hi", "shape")
  p_scale <- stability_plot(stab, u_ref, market, "Adj_Scale",
                            "Adj_Scale_Lo", "Adj_Scale_Hi", "scale")

  p_mef / (p_shape | p_scale)
}

save_figure(build_panel(spx,  "S&P 500"),  "spx_conditional_diagnostics.svg", width = 10, height = 8)
save_figure(build_panel(ftse, "FTSE MIB"), "mib_conditional_diagnostics.svg", width = 10, height = 8)
