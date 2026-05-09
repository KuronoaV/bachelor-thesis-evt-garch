# 12_forecast_mib.R
# Runs the forecasting engine (script 10) on the FTSE MIB over the three
# crisis windows in cfg$crises and saves mib_forecasts.{rds,csv}.

source("scripts/00_setup.R")
source("scripts/10_forecasting_engine.R")

ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

ftse_forecasts <- bind_rows(
  lapply(cfg$crises, function(cr) run_rolling(ftse, "FTSE MIB", cr))
)

saveRDS(ftse_forecasts, file.path(paths$data_proc, "mib_forecasts.rds"))
save_table(ftse_forecasts, "mib_forecasts.csv")
