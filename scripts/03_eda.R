# 03_eda.R
# Full-sample summary statistics, QQ plots vs Normal, and ACF / squared-ACF
# of the daily log-losses (Section 4.3, Figures 4.2 and 4.3).

source("scripts/00_setup.R")

spx  <- readRDS(file.path(paths$data_proc, "spx_log_losses.rds"))
ftse <- readRDS(file.path(paths$data_proc, "ftse_log_losses.rds"))

sample_skewness <- function(x) {
  x <- as.numeric(coredata(x))
  mean((x - mean(x))^3) / sd(x)^3
}
sample_kurtosis <- function(x) {
  x <- as.numeric(coredata(x))
  mean((x - mean(x))^4) / sd(x)^4
}

series <- list(`S&P 500` = spx, `FTSE MIB` = ftse)

full_stats <- tibble(
  Market       = names(series),
  Observations = sapply(series, length),
  Mean         = sapply(series, function(x) mean(coredata(x))),
  SD           = sapply(series, function(x) sd(coredata(x))),
  Skewness     = sapply(series, sample_skewness),
  Kurtosis     = sapply(series, sample_kurtosis)
)
save_table(full_stats, "full_sample_statistics.csv")

qq_panel <- function(x, market) {
  df <- tibble(Loss = as.numeric(coredata(x)))
  ggplot(df, aes(sample = Loss)) +
    stat_qq(colour = palette_market[[market]], alpha = 0.5, shape = 1) +
    stat_qq_line(linetype = "dashed", linewidth = 0.7) +
    labs(title = bquote(.(market) ~ "- QQ vs. Normal"),
         x = "Theoretical Normal quantiles",
         y = "Empirical loss quantiles")
}

# A shared y-limit makes the right-tail asymmetry between SPX and MIB directly
# comparable across panels.
y_lim_qq <- 1.05 * max(abs(c(coredata(spx), coredata(ftse))))
qq_composite <- (qq_panel(spx, "S&P 500") + ylim(-y_lim_qq, y_lim_qq)) |
                (qq_panel(ftse, "FTSE MIB") + ylim(-y_lim_qq, y_lim_qq))
save_figure(qq_composite, "qq_plots.svg", width = 10, height = 5)

acf_frame <- function(x, lag_max = 40) {
  out <- acf(as.numeric(coredata(x)), plot = FALSE, lag.max = lag_max)
  tibble(lag = as.numeric(out$lag[-1]), acf = as.numeric(out$acf[-1]))
}

acf_plot <- function(df, n_obs, title, y_lim) {
  ci <- 1.96 / sqrt(n_obs)
  ggplot(df, aes(lag, acf)) +
    geom_segment(aes(xend = lag, yend = 0), linewidth = 0.3) +
    geom_hline(yintercept = c(-ci, ci), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_x_continuous(breaks = seq(0, max(df$lag), by = 5)) +
    coord_cartesian(ylim = y_lim) +
    labs(title = title, x = "Lag", y = "Autocorrelation")
}

raw_acf <- list(spx = acf_frame(spx),  mib = acf_frame(ftse))
sq_acf  <- list(spx = acf_frame(spx^2), mib = acf_frame(ftse^2))

# Asymmetric y-limit on the squared-loss panels so the slow decay reads clearly
# without wasting space below zero.
raw_lim   <- 0.02 + max(abs(c(raw_acf$spx$acf, raw_acf$mib$acf)))
sq_upper  <- 0.05 + max(c(sq_acf$spx$acf, sq_acf$mib$acf))
sq_lower  <- -sq_upper / 3

acf_composite <- (
  acf_plot(raw_acf$spx, length(spx),  bquote("S&P 500 - ACF of" ~ L[t]),  c(-raw_lim, raw_lim)) |
  acf_plot(raw_acf$mib, length(ftse), bquote("FTSE MIB - ACF of" ~ L[t]), c(-raw_lim, raw_lim))
) / (
  acf_plot(sq_acf$spx,  length(spx),  bquote("S&P 500 - ACF of" ~ L[t]^2),  c(sq_lower, sq_upper)) |
  acf_plot(sq_acf$mib,  length(ftse), bquote("FTSE MIB - ACF of" ~ L[t]^2), c(sq_lower, sq_upper))
)
save_figure(acf_composite, "acf_plots.svg", width = 12, height = 8)
