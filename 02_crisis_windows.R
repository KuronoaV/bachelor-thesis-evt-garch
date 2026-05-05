# 02_crisis_windows.R
#
# Defines the three crisis evaluation windows used in Chapters 4 and 6 and
# produces the visual and tabular evidence that anchors them in the data. The
# canonical dates live in `cfg$crises` (Section 4.2 of the thesis): the
# Dot-com bubble runs from the Nasdaq peak on 24 March 2000 to the S&P 500
# trough on 9 October 2002; the Global Financial Crisis from the S&P peak of
# 9 October 2007 to the trough of 9 March 2009; the COVID-19 crash from the
# 19 February 2020 peak to the 23 March 2020 trough. The same dates feed the
# rolling-forecast scripts (10-12) and the backtest pipeline (13-14), so this
# file is the only place where they are written down.
#
# Two artifacts come out of the script. The first is a per-crisis summary
# table reporting mean, standard deviation, minimum and maximum of the daily
# log-loss inside each window, separately for the two indices. These tables
# read into Section 4.2 and document, in numbers, the asymmetry between the
# slow-burn Dot-com correction, the highly clustered GFC drawdown, and the
# concentrated COVID-19 shock. The second is a 3x2 composite of price
# trajectories with the crisis window shaded; the two-year padding around
# each window lets the reader see the build-up and the recovery, not just the
# eye of the storm.
#
# Numbers are saved as numerics, not as pre-formatted strings: any later
# script that wants to cite "the COVID FTSE drawdown" can read the CSV and
# format it consistently with the rest of Chapter 6 instead of re-parsing
# decimals. The plotting helper takes a single crisis spec and an index name
# and is the only function the reader needs to follow to understand the grid.

source("scripts/00_setup.R")

spx_full  <- readRDS(file.path(paths$data_proc, "spx_full.rds"))  %>% mutate(Market = "S&P 500")
ftse_full <- readRDS(file.path(paths$data_proc, "ftse_full.rds")) %>% mutate(Market = "FTSE MIB")

prices <- bind_rows(spx_full, ftse_full) %>%
  mutate(Market = factor(Market, levels = c("S&P 500", "FTSE MIB")))

summarise_crisis <- function(df, crisis) {
  df %>%
    filter(Date >= crisis$start, Date <= crisis$end) %>%
    group_by(Market) %>%
    summarise(
      Crisis = crisis$label,
      N      = n(),
      Mean   = mean(Log_Loss),
      SD     = sd(Log_Loss),
      Min    = min(Log_Loss),
      Max    = max(Log_Loss),
      .groups = "drop"
    )
}

crisis_stats <- bind_rows(lapply(cfg$crises, summarise_crisis, df = prices))
save_table(crisis_stats, "crisis_summary_statistics.csv")

# y_step lets each panel pick a price grid that suits its own scale; the
# index-specific palette comes from the global setup.
plot_crisis_panel <- function(df, crisis, market, y_step, padding_years = 2) {
  pad        <- padding_years * 365
  plot_start <- crisis$start - pad
  plot_end   <- crisis$end   + pad
  panel      <- df %>% filter(Market == market, Date >= plot_start, Date <= plot_end)

  ggplot(panel, aes(Date, Price)) +
    annotate("rect",
             xmin = crisis$start, xmax = crisis$end,
             ymin = -Inf, ymax = Inf,
             alpha = 0.2, fill = "black") +
    geom_line(colour = palette_market[[market]], linewidth = 0.4) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(breaks = breaks_width(y_step), labels = comma) +
    labs(title = paste(market, "-", crisis$label), x = NULL, y = "Index Price")
}

y_steps <- list(
  Dotcom = list(`S&P 500` = 200, `FTSE MIB` = 5000),
  GFC    = list(`S&P 500` = 200, `FTSE MIB` = 5000),
  Covid  = list(`S&P 500` = 400, `FTSE MIB` = 2500)
)

build_row <- function(key) {
  crisis <- cfg$crises[[key]]
  step   <- y_steps[[key]]
  plot_crisis_panel(prices, crisis, "S&P 500",  step$`S&P 500`) |
    plot_crisis_panel(prices, crisis, "FTSE MIB", step$`FTSE MIB`)
}

composite <- build_row("Dotcom") /
             build_row("GFC")    /
             build_row("Covid")

save_figure(composite, "crisis_price_trajectories.svg", width = 10, height = 12)
