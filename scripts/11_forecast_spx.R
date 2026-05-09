# 11_forecast_spx.R
# Runs the forecasting engine (script 10) on the S&P 500 over the three
# crisis windows in cfg$crises and saves spx_forecasts.{rds,csv}.

source("scripts/00_setup.R")
source("scripts/10_forecasting_engine.R")

spx <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))

spx_forecasts <- bind_rows(
  lapply(cfg$crises, function(cr) run_rolling(spx, "S&P 500", cr))
)

saveRDS(spx_forecasts, file.path(paths$data_proc, "spx_forecasts.rds"))
save_table(spx_forecasts, "spx_forecasts.csv")
