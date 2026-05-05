# 12_forecast_mib.R
#
# Rolling-window driver for the FTSE MIB leg of the Chapter 6 horse race.
# Mirror image of script 11: same forecasting engine (script 10), same
# crisis windows (cfg$crises), same per-day output schema, only the input
# series and the output filenames change. Keeping the two drivers as thin
# wrappers around a single shared engine is what guarantees that any
# difference between the SPX and MIB columns in Chapter 6 is a difference
# in the data, not a difference in the code path.
#
# The one place where SPX and MIB do diverge is the available history at
# the start of the Dot-com window. The FTSE MIB series begins in January
# 1992 (with the MIB 30 spliced into FTSE MIB), which gives roughly 2050
# trading days by 1 March 2000 — short of the cfg$windows$uevt = 2500
# target. The dynamic Unconditional EVT logic baked into run_rolling
# (script 10) handles this transparently: it falls back to whatever
# history is available and emits a warning whenever the window drops
# below cfg$windows$uevt_min = 1500. Those warnings are informative and
# should be reviewed when the script is rerun, but they are expected for
# the FTSE MIB and do not indicate a failure.

source("scripts/00_setup.R")
source("scripts/10_forecasting_engine.R")

ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

ftse_forecasts <- bind_rows(
  lapply(cfg$crises, function(cr) run_rolling(ftse, "FTSE MIB", cr))
)

saveRDS(ftse_forecasts, file.path(paths$data_proc, "mib_forecasts.rds"))
save_table(ftse_forecasts, "mib_forecasts.csv")
