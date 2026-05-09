# 15_rolling_shape_parameter.R
# Day-by-day GPD shape estimate xi over the GFC window for the S&P 500,
# computed in parallel for the Unconditional EVT model (raw losses) and the
# Conditional EVT model (standardised residuals).

source("scripts/00_setup.R")

if (!requireNamespace("rugarch",  quietly = TRUE)) install.packages("rugarch")
if (!requireNamespace("extRemes", quietly = TRUE)) install.packages("extRemes")
suppressPackageStartupMessages({
  library(rugarch)
  library(extRemes)
})

spx <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))

crisis     <- cfg$crises$GFC
loss_vec   <- as.numeric(coredata(spx))
date_vec   <- index(spx)
eval_dates <- date_vec[date_vec >= crisis$start & date_vec <= crisis$end]

garch_spec <- ugarchspec(
  variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
  distribution.model = "norm"
)

shape_estimate <- function(x, u_quantile = cfg$u_quantile) {
  u   <- quantile(x, probs = u_quantile, names = FALSE)
  fit <- tryCatch(
    suppressWarnings(fevd(x, threshold = u, type = "GP", method = "MLE")),
    error = function(e) NULL
  )
  if (is.null(fit)) NA_real_ else fit$results$par[["shape"]]
}

shape_uncond <- function(losses) shape_estimate(losses)

shape_cond <- function(losses) {
  fit <- tryCatch(
    suppressWarnings(ugarchfit(spec = garch_spec, data = -losses, solver = "hybrid")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  z_loss <- -as.numeric(residuals(fit, standardize = TRUE))
  shape_estimate(z_loss)
}

cat(sprintf("[GFC | %s] %d days\n", "S&P 500", length(eval_dates)))

rows <- lapply(seq_along(eval_dates), function(j) {
  t <- which(date_vec == eval_dates[j])
  w_uevt <- loss_vec[(t - cfg$windows$uevt):(t - 1L)]
  w_cevt <- loss_vec[(t - cfg$windows$cevt):(t - 1L)]

  if (j %% 50L == 0L) cat(sprintf("  %d / %d\n", j, length(eval_dates)))

  tibble(
    Date             = date_vec[t],
    Xi_Unconditional = shape_uncond(w_uevt),
    Xi_Conditional   = shape_cond(w_cevt)
  )
})

shape_df <- bind_rows(rows)
save_table(shape_df, "rolling_shape_gfc.csv")

plot_df <- shape_df %>%
  pivot_longer(cols = c(Xi_Unconditional, Xi_Conditional),
               names_to = "Model", values_to = "Xi") %>%
  mutate(Model = factor(
    Model,
    levels = c("Xi_Unconditional", "Xi_Conditional"),
    labels = c("Unconditional EVT (n=2500, raw losses)",
               "Conditional EVT (n=1000, residuals)")
  ))

p <- ggplot(plot_df, aes(Date, Xi, colour = Model)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(
    "Unconditional EVT (n=2500, raw losses)" = palette_model[["Unconditional EVT"]],
    "Conditional EVT (n=1000, residuals)"    = palette_model[["Conditional EVT"]]
  )) +
  labs(title = bquote("Rolling GPD shape" ~ xi ~ "over the GFC (S&P 500)"),
       x = NULL, y = expression(xi), colour = NULL)

save_figure(p, "rolling_shape_comparison.svg", width = 10, height = 5)
