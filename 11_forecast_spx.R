# 11_forecast_spx.R
#
# Rolling-window driver for the S&P 500 leg of the Chapter 6 horse race.
# Sources the forecasting engine of script 09 and walks each of the three
# crisis windows defined in cfg$crises (Section 4.2: Dot-com, GFC, COVID),
# producing one row per evaluation day with the realised log-loss and the
# eight model forecasts (VaR and ES at the FRTB levels for each of the four
# competing models). The output is the input for the backtests of script 12
# and the result tables and figures of script 13.
#
# Two practical notes. First, the S&P 500 history begins in January 1990, so
# by the start of the Dot-com window in 2000 the series carries the full
# 2500-day history that the Unconditional EVT model targets; no dynamic
# fallback fires for SPX, and any UEVT short-history warnings should be
# read as an unintended config drift rather than an expected event.
# Second, the evaluation window is padded by 30 calendar days on each side
# of the canonical crisis dates by the rolling driver — see the header of
# script 09 for why — and that padding is what gives the COVID episode
# enough days to support the Kupiec, Christoffersen and Acerbi-Szekely
# tests of Chapter 6 with non-degenerate sample sizes.

source("scripts/00_setup.R")
source("scripts/10_forecasting_engine.R")

spx <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))

spx_forecasts <- bind_rows(
  lapply(cfg$crises, function(cr) run_rolling(spx, "S&P 500", cr))
)

saveRDS(spx_forecasts, file.path(paths$data_proc, "spx_forecasts.rds"))
save_table(spx_forecasts, "spx_forecasts.csv")
