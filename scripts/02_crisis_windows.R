# 02_crisis_windows.R
# Per-crisis summary statistics (Table 4.2) and the price-trajectory composite
# with each crisis window shaded (Figure 4.1).

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
