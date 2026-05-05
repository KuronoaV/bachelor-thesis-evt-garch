# 10_forecasting_engine.R
#
# This file is the forecasting engine of the thesis: it defines the four
# competing one-day-ahead VaR/ES forecasters of Section 3.6, plus the rolling
# driver that script 11 (S&P 500) and script 12 (FTSE MIB) call to produce
# the Chapter 6 evaluation set. Nothing executes at the top level — the file
# is a library, sourced by the per-market drivers and never run on its own.
# Putting the forecasters and the driver in one place removes the duplication
# the original code carried between SPX and MIB and ensures every forecast in
# the thesis comes from the same code path.
#
# The four forecasters share the same signature: each takes a numeric vector
# of past losses and returns a named vector c(VaR, ES) at the FRTB levels
# (alpha_var = 0.99, alpha_es = 0.975) defined in the cfg list. Gaussian
# Parametric assumes Normality and reads VaR off qnorm and ES off the
# standard hazard-ratio formula. Historical Simulation reads VaR off the
# empirical 99% quantile and ES off the mean of losses strictly above the
# 97.5% quantile. Unconditional EVT fits a GPD by MLE above the 90th
# percentile and applies the Pickands-Balkema-de Haan tail-quantile formula.
# Conditional EVT runs the McNeil-Frey two-step procedure: an
# AR(1)-GJR-GARCH(1,1) filter followed by a GPD on the standardised
# innovations, with the next-day forecast assembled as
# mu_loss(t+1) + sigma(t+1) * (innovation-tail quantile).
#
# Three implementation choices are worth flagging here. First, the GARCH
# step is fitted on returns r_t = -L_t to keep the leverage parameter
# positive, in line with Section 3.5; the residuals are flipped back to
# loss form before the GPD step so the upper tail is the loss tail. Second,
# the GPD VaR/ES helper handles the two numerical edge cases that bite the
# textbook formula: it uses the limiting exponential form when |shape| is
# below 1e-6 and returns NA for ES when shape >= 1, since the tail mean is
# then infinite. Third, every fit is wrapped in a tryCatch that returns NA
# rather than throwing, because the rolling driver hits ~150 windows per
# crisis and one optimiser failure would otherwise abort the whole run.
#
# The rolling driver iterates a single crisis spec (label, slug, start, end)
# and produces one row per evaluation day with the realised loss and the
# eight model forecasts. The evaluation window is padded by `pad_days`
# calendar days on each side of the crisis dates so the backtests of
# Chapter 6 get enough sample to be powered, especially the COVID episode
# whose canonical window is only a month long. The Unconditional EVT
# rolling window is dynamic: it targets cfg$windows$uevt = 2500 days but
# falls back to whatever history is available, with a warning if it drops
# below cfg$windows$uevt_min = 1500. This is what lets the same driver
# serve both the long-history SPX series and the shorter FTSE MIB series.

source("scripts/00_setup.R")

if (!requireNamespace("rugarch",  quietly = TRUE)) install.packages("rugarch")
if (!requireNamespace("extRemes", quietly = TRUE)) install.packages("extRemes")
suppressPackageStartupMessages({
  library(rugarch)
  library(extRemes)
})

# Numerical companion to the four forecasters: the Pickands-Balkema-de Haan
# tail-quantile formula with the shape ~ 0 limit handled cleanly and ES
# returned as NA when shape >= 1.
gpd_var_es <- function(scale, shape, threshold, n_total, n_exc, alpha_var, alpha_es) {
  if (is.na(scale) || is.na(shape) || scale <= 0)
    return(c(VaR = NA_real_, ES = NA_real_))

  pot_var <- function(alpha) {
    q <- (n_total / n_exc) * (1 - alpha)
    if (abs(shape) < 1e-6) threshold + scale * (-log(q))
    else                   threshold + (scale / shape) * (q^(-shape) - 1)
  }

  v_var <- pot_var(alpha_var)
  if (shape >= 1) {
    es <- NA_real_
  } else {
    v_es <- pot_var(alpha_es)
    es <- if (abs(shape) < 1e-6) v_es + scale
          else                   (v_es + scale - shape * threshold) / (1 - shape)
  }

  c(VaR = unname(v_var), ES = unname(es))
}

forecast_gaussian <- function(losses,
                              alpha_var = cfg$alpha_var,
                              alpha_es  = cfg$alpha_es) {
  mu    <- mean(losses)
  sigma <- sd(losses)
  c(VaR = mu + sigma * qnorm(alpha_var),
    ES  = mu + sigma * dnorm(qnorm(alpha_es)) / (1 - alpha_es))
}

forecast_hs <- function(losses,
                        alpha_var = cfg$alpha_var,
                        alpha_es  = cfg$alpha_es) {
  v   <- as.numeric(quantile(losses, probs = alpha_var, names = FALSE))
  thr <- as.numeric(quantile(losses, probs = alpha_es,  names = FALSE))
  c(VaR = v, ES = mean(losses[losses > thr]))
}

forecast_uncond_evt <- function(losses,
                                alpha_var  = cfg$alpha_var,
                                alpha_es   = cfg$alpha_es,
                                u_quantile = cfg$u_quantile,
                                min_exc    = 10) {
  n  <- length(losses)
  u  <- quantile(losses, probs = u_quantile, names = FALSE)
  Nu <- sum(losses > u)
  if (Nu < min_exc) return(c(VaR = NA_real_, ES = NA_real_))

  fit <- tryCatch(
    suppressWarnings(fevd(losses, threshold = u, type = "GP", method = "MLE")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(c(VaR = NA_real_, ES = NA_real_))

  par <- fit$results$par
  gpd_var_es(par[["scale"]], par[["shape"]], u, n, Nu, alpha_var, alpha_es)
}

forecast_cond_evt <- function(losses,
                              alpha_var  = cfg$alpha_var,
                              alpha_es   = cfg$alpha_es,
                              u_quantile = cfg$u_quantile,
                              min_exc    = 10) {

  spec <- ugarchspec(
    variance.model     = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model         = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = "norm"
  )

  fit_g <- tryCatch(
    suppressWarnings(ugarchfit(spec = spec, data = -losses, solver = "hybrid")),
    error = function(e) NULL
  )
  if (is.null(fit_g)) return(c(VaR = NA_real_, ES = NA_real_))

  fcst <- tryCatch(ugarchforecast(fit_g, n.ahead = 1), error = function(e) NULL)
  if (is.null(fcst)) return(c(VaR = NA_real_, ES = NA_real_))

  mu_r <- as.numeric(fitted(fcst))
  s_n  <- as.numeric(sigma(fcst))
  if (any(is.na(c(mu_r, s_n))) || s_n <= 0)
    return(c(VaR = NA_real_, ES = NA_real_))

  z_loss <- -as.numeric(residuals(fit_g, standardize = TRUE))

  n  <- length(z_loss)
  u  <- quantile(z_loss, probs = u_quantile, names = FALSE)
  Nu <- sum(z_loss > u)
  if (Nu < min_exc) return(c(VaR = NA_real_, ES = NA_real_))

  fit_e <- tryCatch(
    suppressWarnings(fevd(z_loss, threshold = u, type = "GP", method = "MLE")),
    error = function(e) NULL
  )
  if (is.null(fit_e)) return(c(VaR = NA_real_, ES = NA_real_))

  par <- fit_e$results$par
  inn <- gpd_var_es(par[["scale"]], par[["shape"]], u, n, Nu, alpha_var, alpha_es)

  mu_loss_next <- -mu_r
  c(VaR = unname(mu_loss_next + s_n * inn[["VaR"]]),
    ES  = unname(mu_loss_next + s_n * inn[["ES"]]))
}

# Walks day by day through the (padded) crisis window and produces one
# tibble row per evaluation day. UEVT degrades from cfg$windows$uevt to
# whatever history is available and warns if it falls below uevt_min.
run_rolling <- function(losses, market, crisis,
                        pad_days = 0,
                        windows  = cfg$windows) {
  loss_vec <- as.numeric(coredata(losses))
  date_vec <- index(losses)

  start_date <- crisis$start - pad_days
  end_date   <- crisis$end   + pad_days
  eval_idx   <- which(date_vec >= start_date & date_vec <= end_date)

  min_history <- max(windows$gauss_hs, windows$cevt)
  eval_idx    <- eval_idx[eval_idx > min_history]

  if (length(eval_idx) == 0L) {
    warning(sprintf("%s [%s]: no eval days with sufficient history",
                    market, crisis$label))
    return(tibble())
  }

  cat(sprintf("[%s | %s] %d days\n", market, crisis$label, length(eval_idx)))

  rows <- lapply(seq_along(eval_idx), function(j) {
    t   <- eval_idx[j]
    avl <- t - 1L

    w_short <- loss_vec[(t - windows$gauss_hs):(t - 1L)]
    w_cevt  <- loss_vec[(t - windows$cevt):(t - 1L)]

    uevt_n <- min(windows$uevt, avl)
    if (uevt_n < windows$uevt_min) {
      warning(sprintf("UEVT history short on %s: %d days (<%d)",
                      date_vec[t], uevt_n, windows$uevt_min))
    }
    w_uevt <- loss_vec[(t - uevt_n):(t - 1L)]

    g  <- forecast_gaussian(w_short)
    h  <- forecast_hs(w_short)
    u  <- forecast_uncond_evt(w_uevt)
    ce <- forecast_cond_evt(w_cevt)

    if (j %% 50L == 0L) cat(sprintf("  %d / %d\n", j, length(eval_idx)))

    tibble(
      Market        = market,
      Crisis        = crisis$label,
      Date          = date_vec[t],
      Realized_Loss = loss_vec[t],
      Gauss_VaR     = g[["VaR"]],  Gauss_ES = g[["ES"]],
      HS_VaR        = h[["VaR"]],  HS_ES    = h[["ES"]],
      UEVT_VaR      = u[["VaR"]],  UEVT_ES  = u[["ES"]],
      CEVT_VaR      = ce[["VaR"]], CEVT_ES  = ce[["ES"]]
    )
  })

  bind_rows(rows)
}
