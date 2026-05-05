# Extreme Value Theory and Tail Risk Forecasting

This repository contains the complete empirical R codebase for my Master's thesis. The project evaluates the forecasting performance of Value-at-Risk (VaR) and Expected Shortfall (ES) models during severe financial market crises (e.g., the Dot-com bubble, the Global Financial Crisis, and the COVID-19 crash). 

Specifically, this codebase implements and backtests the two-step **Conditional Extreme Value Theory (EVT)** methodology proposed by McNeil and Frey (2000), utilizing an AR(1)-GJR-GARCH(1,1) filter to model volatility clustering, followed by the application of the Peaks-Over-Threshold (POT) method on the standardized innovations.

## ⚠️ Data Availability Statement

The empirical analysis in this thesis utilizes proprietary daily index data (S&P 500 and FTSE MIB) sourced from **Bloomberg**. Under the Bloomberg Terminal terms of service and strict data redistribution agreements, the raw and processed datasets cannot be publicly shared or hosted on GitHub. 

Therefore, the `data/` directory (containing the `.rds` and `.csv` files) has been intentionally excluded from this repository to ensure full legal compliance. The R scripts provided here demonstrate the complete, reproducible econometric pipeline and can be executed by substituting any appropriately formatted financial time-series data.

## Repository Structure

The codebase is highly modular, separating data preparation, in-sample estimation, out-of-sample forecasting, and statistical backtesting into distinct operational scripts:

* **`00_setup.R`**: Global configurations, library management, and plot aesthetics.
* **`01` to `03`**: Exploratory data analysis, extraction of crisis windows, and generation of baseline stylized facts.
* **`04` to `08a`**: In-sample Chapter 5 diagnostics, including unconditional GPD fitting, GARCH residual testing, and mean-excess function stability checks.
* **`10_forecasting_engine.R`**: The core mathematical engine containing the daily rolling algorithms for Gaussian Parametric, Historical Simulation, Unconditional EVT, and Conditional EVT models.
* **`11_forecast_spx.R` & `12_forecast_mib.R`**: The out-of-sample driver scripts that execute the rolling forecasts over the designated crisis windows.
* **`13_backtesting_engine.R`**: A rigorous statistical backtesting library implementing Kupiec Unconditional Coverage, Christoffersen Independence/Conditional Coverage, Acerbi-Szekely $Z_2$ tests, and a 10,000-path Monte Carlo parametric bootstrap for finite-sample inference.
* **`14_chapter6_results.R`**: The orchestrator that ingests the forecasts, runs the backtests, generates the final evaluation matrices, and plots the McNeil-style VaR/ES violation grids.

## Methodological Notes

* **No-Look-Ahead Bias:** The rolling forecasting engine strictly slices data up to $t-1$ when projecting risk metrics for day $t$.
* **Finite-Sample Robustness:** Because crisis windows are notoriously short (e.g., ~22 trading days for the COVID-19 crash), standard asymptotic $\chi^2$ backtests lose statistical power. This pipeline utilizes a Monte Carlo bootstrap engine (`13_backtesting_engine.R`) to compute exact finite-sample $p$-values.
