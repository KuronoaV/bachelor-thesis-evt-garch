# 04_pot_diagnostics.R
# Mean-excess function and threshold-stability plots for the GPD shape and
# modified scale, on the raw losses (Section 4.4, Figures 4.4 and 4.5).

source("scripts/00_setup.R")

if (!requireNamespace("extRemes", quietly = TRUE)) install.packages("extRemes")
suppressPackageStartupMessages(library(extRemes))

spx  <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))
ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

mef_grid <- function(x, qs = seq(0.50, 0.995, by = 0.005)) {
  x  <- as.numeric(coredata(x))
  us <- quantile(x, probs = qs, names = FALSE)
  tibble(
    Quantile  = qs,
    Threshold = us,
    MEF       = sapply(us, function(u) mean(x[x > u] - u))
  )
}

# Wald CIs straight off the Hessian. Drops thresholds with too few
# exceedances or where the optimiser fails to converge. The CI for the
# modified scale beta* = beta - xi*u uses the delta method so that the
# covariance between beta and xi is propagated correctly:
#   Var(beta*) = Var(beta) + u^2 Var(xi) - 2u Cov(beta, xi)
stability_grid <- function(x, qs = seq(0.50, 0.985, by = 0.005), min_exc = 25) {
  x  <- as.numeric(coredata(x))
  us <- quantile(x, probs = qs, names = FALSE)

  rows <- lapply(seq_along(us), function(i) {
    u <- us[i]
    if (sum(x > u) < min_exc) return(NULL)
    fit <- tryCatch(
      suppressWarnings(fevd(x, threshold = u, type = "GP", method = "MLE")),
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
      Quantile      = qs[i],
      Threshold     = u,
      Shape         = shape,
      Shape_Lo      = shape - 1.96 * shape_se,
      Shape_Hi      = shape + 1.96 * shape_se,
      Adj_Scale     = beta_star,
      Adj_Scale_Lo  = beta_star - 1.96 * beta_star_se,
      Adj_Scale_Hi  = beta_star + 1.96 * beta_star_se
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
    labs(title = bquote(.(market) ~ "- Mean excess function"),
         x = "Threshold (u)", y = "Mean excess")
}

stability_plot <- function(df, u_ref, market, value, lo, hi, kind) {
  if (kind == "shape") {
    title_expr <- bquote(.(market) ~ "-" ~ xi ~ "stability")
    y_expr     <- expression(xi)
  } else {
    title_expr <- bquote(.(market) ~ "-" ~ beta * "*" ~ "stability")
    y_expr     <- expression(beta * "*")
  }
  ggplot(df, aes(Threshold, .data[[value]])) +
    geom_ribbon(aes(ymin = .data[[lo]], ymax = .data[[hi]]),
                fill = palette_market[[market]], alpha = 0.15) +
    geom_line(linewidth = 0.5, colour = palette_market[[market]]) +
    geom_vline(xintercept = u_ref, linetype = "dashed", colour = "red") +
    labs(title = title_expr, x = "Threshold (u)", y = y_expr)
}

build_panel <- function(x, market) {
  u_ref <- quantile(as.numeric(coredata(x)), probs = cfg$u_quantile, names = FALSE)
  mef   <- mef_grid(x)
  stab  <- stability_grid(x)

  p_mef   <- mef_plot(mef, u_ref, market)
  p_shape <- stability_plot(stab, u_ref, market, "Shape",
                            "Shape_Lo", "Shape_Hi", "shape")
  p_scale <- stability_plot(stab, u_ref, market, "Adj_Scale",
                            "Adj_Scale_Lo", "Adj_Scale_Hi", "scale")

  p_mef / (p_shape | p_scale)
}

save_figure(build_panel(spx,  "S&P 500"),  "spx_pot_diagnostics.svg", width = 10, height = 8)
save_figure(build_panel(ftse, "FTSE MIB"), "mib_pot_diagnostics.svg", width = 10, height = 8)
