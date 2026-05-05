# 15_rolling_shape_parameter.R
#
# Day-by-day GPD shape estimate xi over the Global Financial Crisis
# evaluation window, computed in parallel for the Unconditional EVT model
# (raw losses, n = cfg$windows$uevt) and the Conditional EVT model
# (standardised innovations, n = cfg$windows$cevt). The point of the
# figure is the contrast: the unconditional shape drifts noticeably as
# the GFC enters the rolling window — the heaviest tail observations get
# included, then eventually exit out the back of the window — while the
# conditional shape tracks the residuals after the GARCH filter has
# already absorbed the volatility component, and is comparatively stable.
# This is the empirical anchor for the discussion in Section 7.3 of the
# stability advantages that the McNeil-Frey two-step procedure buys over
# its unconditional cousin.
#
# The GFC window is the natural laboratory for this comparison because it
# is the longest of the three crises, so the rolling drift in the
# unconditional estimate has time to develop and decay; the Dot-com
# window also covers the right span but the contrast is muted on the
# S&P 500, and the COVID window is too short to be visually informative.
# The script is therefore deliberately scoped to one market (S&P 500) and
# one crisis (GFC); reproducing the same exercise for the FTSE MIB or
# for a different crisis is a one-line change to the inputs at the top.
#
# Implementation choices follow the rest of the pipeline. Window sizes
# come from cfg$windows; the GFC dates come from cfg$crises$GFC; the
# threshold rule is cfg$u_quantile applied separately to the raw losses
# and to the residual innovations on each day. The AR(1)-GJR-GARCH(1,1)
# is fitted on returns and the residuals are flipped back to loss form,
# matching the convention introduced in script 06. Each daily fit is
# wrapped in a tryCatch that returns NA for that day's xi if the
# optimiser refuses to converge, since one bad day should not abort the
# loop over hundreds of evaluation points.

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
