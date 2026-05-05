# 14_chapter6_results.R
#
# Orchestrator for the Chapter 6 results: takes the per-day forecasts produced
# by scripts 11 (S&P 500) and 12 (FTSE MIB), feeds them through the backtest
# library of script 13, and writes the comprehensive backtest table that
# Chapter 6 reads from. It also produces the per-crisis × per-market McNeil
# grid plots that visualise where each model's VaR and ES forecasts land
# relative to the realised losses, with the violation taxonomy that Section
# 6.1 leans on for the qualitative discussion. Nothing modelling-related is
# defined in this script — it is purely a glue layer that executes the
# library code and saves the artifacts.
#
# The backtest table is built one Market × Crisis × Model row at a time and
# carries every statistic that Chapter 6 cites: expected and observed
# violation counts, Basel zone, both Kupiec p-values (asymptotic and Monte
# Carlo), the Christoffersen LR_ind / LR_cc pair with their p-values, the
# Acerbi-Szekely Z2 and its zone. Numerics are saved at full precision and
# never pre-rounded, in line with the convention of the rest of the thesis
# pipeline: the LaTeX-side formatters are responsible for choosing how many
# decimals appear in the printed tables, not this script. A second slim CSV
# (traffic_light_summary.csv) projects only the two zone classifications,
# which is the at-a-glance view that Section 6.4 uses as a wrap-up table.
#
# The McNeil grids show 99% VaR and 97.5% ES as two grey lines per model,
# the realised log-losses as vertical spikes from the baseline, and the
# violations as coloured points: blue squares for VaR-only breaches and
# orange triangles for the more severe breaches that pierce the ES line as
# well. Six figures come out of the loop — three crises times two markets —
# and the file naming follows the cfg$crises slug so that the LaTeX side
# can reference them mechanically. The four model panels share a common
# y-axis within each figure so that Chapter 6 can compare model
# conservatism without having to read off different scales.

source("scripts/00_setup.R")
source("scripts/13_backtesting_engine.R")

if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
suppressPackageStartupMessages(library(tidyr))

spx_forecasts <- readRDS(file.path(paths$data_proc, "spx_forecasts.rds"))
mib_forecasts <- readRDS(file.path(paths$data_proc, "mib_forecasts.rds"))
all_forecasts <- bind_rows(spx_forecasts, mib_forecasts)

models <- c("Gauss", "HS", "UEVT", "CEVT")

backtest_one <- function(df_sub, alpha_var = cfg$alpha_var, alpha_es = cfg$alpha_es) {
  L <- df_sub$Realized_Loss
  T <- nrow(df_sub)
  
  bind_rows(lapply(models, function(mod) {
    var_f <- df_sub[[paste0(mod, "_VaR")]]
    es_f  <- df_sub[[paste0(mod, "_ES")]]
    
    kup    <- kupiec_test(L, var_f, alpha_var)
    chris  <- christoffersen_test(L, var_f, alpha_var)
    mc     <- var_backtest_mc(L, var_f, alpha_var)
    z2     <- acerbi_szekely_z2(L, var_f, es_f, alpha_es)
    basel  <- basel_traffic_light(T, kup$violations, alpha_var)
    
    tibble(
      Market              = unique(df_sub$Market),
      Crisis              = unique(df_sub$Crisis),
      Model               = mod,
      N_Days              = T,
      Expected_Violations = kup$expected,
      Observed_Violations = kup$violations,
      Basel_Zone          = basel,
      Kupiec_LRuc         = kup$LRuc,
      Kupiec_p            = kup$p_value,
      Kupiec_MC_p         = mc$MC_p_uc,
      Christo_LRind       = chris$LRind,
      Christo_p_ind       = chris$p_value_ind,
      Christo_MC_p_ind    = mc$MC_p_ind,
      Christo_LRcc        = chris$LRcc,
      Christo_p_cc        = chris$p_value_cc,
      Christo_MC_p_cc     = mc$MC_p_cc,
      Acerbi_Z2           = z2$Z2,
      Acerbi_Zone         = z2$Classification
    )
  }))
}

run_backtests <- function(df) {
  df %>%
    group_by(Market, Crisis) %>%
    group_split() %>%
    lapply(backtest_one) %>%
    bind_rows()
}

backtest_results <- run_backtests(all_forecasts)
save_table(backtest_results, "chapter6_backtest_results.csv")

save_table(
  backtest_results %>% dplyr::select(Market, Crisis, Model, Basel_Zone, Acerbi_Zone),
  "traffic_light_summary.csv"
)

# Define clean labels for the facet headers
model_labels <- c(
  "Gauss" = "Gaussian Parametric",
  "HS" = "Historical Simulation",
  "UEVT" = "Unconditional EVT",
  "CEVT" = "Conditional EVT"
)

# McNeil-style grid: VaR and ES as lines, realised losses as spikes,
# violations as coloured markers split by severity.
mcneil_grid <- function(df, market, crisis) {
  sub <- df %>% filter(Market == market, Crisis == crisis$label)
  
  long <- sub %>%
    dplyr::select(Date, Realized_Loss,
                  Gauss_VaR, Gauss_ES, HS_VaR, HS_ES,
                  UEVT_VaR, UEVT_ES, CEVT_VaR, CEVT_ES) %>%
    pivot_longer(cols = -c(Date, Realized_Loss),
                 names_to = c("Model", ".value"),
                 names_sep = "_") %>%
    mutate(
      Violation = case_when(
        Realized_Loss > ES                                          ~ "VaR + ES",
        Realized_Loss > VaR & Realized_Loss <= ES                   ~ "VaR only",
        TRUE                                                        ~ NA_character_
      ),
      Violation_y = ifelse(!is.na(Violation), Realized_Loss, NA_real_),
      Model       = factor(model_labels[Model], levels = model_labels)
    )
  
  ggplot(long, aes(Date)) +
    geom_linerange(aes(ymin = 0, ymax = Realized_Loss),
                   colour = "grey70", linewidth = 0.4) +
    geom_line(aes(y = VaR, linetype = "99% VaR"),  colour = "grey20", linewidth = 0.6) +
    geom_line(aes(y = ES,  linetype = "97.5% ES"), colour = "grey20", linewidth = 0.6) +
    geom_point(aes(y = Violation_y, shape = Violation, colour = Violation),
               size = 1.4, stroke = 0.8, na.rm = TRUE) +
    scale_linetype_manual(name = "Forecasts",
                          values = c("99% VaR" = "11", "97.5% ES" = "31")) +
    scale_shape_manual(name = "Violations",
                       values = c("VaR + ES" = 2, "VaR only" = 0),
                       na.translate = FALSE) +
    scale_colour_manual(name = "Violations",
                        values = c("VaR + ES" = "#E66101", "VaR only" = "#1f77b4"),
                        na.translate = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)),
                       breaks = pretty_breaks(n = 8)) +
    facet_wrap(~ Model, ncol = 2) +
    labs(title = bquote(.(market) ~ "-" ~ .(crisis$label)), x = NULL, y = "Loss")
}

market_slug <- c("S&P 500" = "spx", "FTSE MIB" = "mib")

for (m in names(market_slug)) {
  for (k in names(cfg$crises)) {
    cr <- cfg$crises[[k]]
    p  <- mcneil_grid(all_forecasts, m, cr)
    save_figure(p,
                sprintf("%s_%s_forecasts.svg", market_slug[[m]], cr$slug),
                width = 12, height = 10)
  }
}