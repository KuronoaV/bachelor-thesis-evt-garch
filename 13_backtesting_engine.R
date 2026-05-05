# 13_backtesting_engine.R
#
# Backtesting library for the rolling forecasts produced by scripts 11 and
# 12. Like script 10 (the forecasting engine), this file is sourced by other
# scripts and never run on its own; it executes nothing at the top level
# beyond loading dependencies. Five test families are implemented and follow
# the order in which they enter the discussion of Chapter 6: the Kupiec
# unconditional-coverage likelihood-ratio test for the violation rate
# (Section 4.5.1), the Christoffersen independence and conditional-coverage
# tests for clustering (Section 4.5.2), a Monte Carlo bootstrap that
# replaces the asymptotic chi-square reference distribution by a simulated
# small-sample one (Section 4.5.3), the Acerbi-Szekely Z2 statistic for
# Expected Shortfall (Section 4.5.4), and the dynamic Basel traffic-light
# classification scaled to the actual sample size by the Binomial CDF
# (Section 4.5.5).
#
# The point of having both the asymptotic Kupiec p-value and the Monte
# Carlo p-value is that crisis windows are short — the COVID episode runs
# for fewer than fifty trading days even after padding — and the chi-square
# null distribution is a poor approximation in that regime. The MC version
# keeps the same likelihood-ratio statistic but replaces the reference
# distribution by 10,000 simulated Bernoulli sequences of the same length
# under the null. The Z2 statistic and the Basel scaling complete the
# picture for ES and for the regulatory-style visual triage; the latter is
# read as a colour, not as a p-value, so its role in the thesis is more
# illustrative than inferential.
#
# Implementation notes. The boundary cases of the LR_uc statistic
# (zero violations or all-violations) are handled by direct algebra rather
# than by the generic log-likelihood difference, since the latter divides
# by zero. The Christoffersen log-likelihoods are evaluated with a small
# helper xlogp(x, p) := x * log(p) on the cells where x > 0 and 0
# otherwise; this implements the standard 0 * log(0) := 0 convention
# without numerical fudging. Reproducibility of the bootstrap is inherited
# from the global set.seed(cfg$mc$seed) call in 00_setup.R; do not re-seed
# inside this file or the per-run results will drift if the call order
# changes upstream.

source("scripts/00_setup.R")

xlogp <- function(x, p) if (x == 0L) 0 else x * log(p)

lruc_stat <- function(violations, t_total, p) {
  if (t_total == 0L) return(NA_real_)
  if (violations == 0L)        return(-2 * t_total * log(1 - p))
  if (violations == t_total)   return(-2 * t_total * log(p))
  pi_hat <- violations / t_total
  log_L0 <- violations * log(p) + (t_total - violations) * log(1 - p)
  log_LA <- violations * log(pi_hat) + (t_total - violations) * log(1 - pi_hat)
  -2 * (log_L0 - log_LA)
}

kupiec_test <- function(realized_losses, var_forecasts, alpha = cfg$alpha_var) {
  p   <- 1 - alpha
  vio <- realized_losses > var_forecasts
  vio <- vio[!is.na(vio)]
  T   <- length(vio)
  x   <- sum(vio)
  
  if (T == 0L) return(list(LRuc = NA_real_, p_value = NA_real_,
                           violations = 0L, expected = 0))
  
  LRuc <- lruc_stat(x, T, p)
  list(LRuc       = LRuc,
       p_value    = 1 - pchisq(LRuc, df = 1),
       violations = x,
       expected   = p * T)
}

christoffersen_test <- function(realized_losses, var_forecasts,
                                alpha = cfg$alpha_var) {
  vio <- as.numeric(realized_losses > var_forecasts)
  vio <- vio[!is.na(vio)]
  T   <- length(vio)
  
  if (T < 2L) {
    return(list(LRind = NA_real_, p_value_ind = NA_real_,
                LRcc  = NA_real_, p_value_cc  = NA_real_))
  }
  
  v_today    <- vio[-T]
  v_tomorrow <- vio[-1]
  
  n00 <- sum(v_today == 0 & v_tomorrow == 0)
  n01 <- sum(v_today == 0 & v_tomorrow == 1)
  n10 <- sum(v_today == 1 & v_tomorrow == 0)
  n11 <- sum(v_today == 1 & v_tomorrow == 1)
  
  # If either row is empty we cannot estimate a transition probability.
  if ((n00 + n01) == 0L || (n10 + n11) == 0L) {
    return(list(LRind = NA_real_, p_value_ind = NA_real_,
                LRcc  = NA_real_, p_value_cc  = NA_real_))
  }
  
  pi01 <- n01 / (n00 + n01)
  pi11 <- n11 / (n10 + n11)
  pi2  <- (n01 + n11) / (T - 1)
  
  log_L0 <- xlogp(n00 + n10, 1 - pi2) + xlogp(n01 + n11, pi2)
  log_LA <- xlogp(n00, 1 - pi01) + xlogp(n01, pi01) +
    xlogp(n10, 1 - pi11) + xlogp(n11, pi11)
  
  LRind <- max(0, -2 * (log_L0 - log_LA))
  
  kup  <- kupiec_test(realized_losses, var_forecasts, alpha)
  LRcc <- kup$LRuc + LRind
  
  list(LRind       = LRind,
       p_value_ind = 1 - pchisq(LRind, df = 1),
       LRcc        = LRcc,
       p_value_cc  = 1 - pchisq(LRcc, df = 2))
}

# Replaces the asymptotic chi-square reference with a simulated null
# for the Kupiec, Christoffersen Independence, and Conditional Coverage tests.
var_backtest_mc <- function(realized_losses, var_forecasts, 
                            alpha = cfg$alpha_var, 
                            n_sim = cfg$mc$n_sim) {
  
  vio <- as.numeric(realized_losses > var_forecasts)
  vio <- vio[!is.na(vio)]
  t_total <- length(vio)
  
  if (t_total < 2L) {
    return(list(MC_p_uc = NA_real_, MC_p_ind = NA_real_, MC_p_cc = NA_real_))
  }
  
  p <- 1 - alpha
  
  obs_kupiec <- kupiec_test(realized_losses, var_forecasts, alpha)
  obs_christ <- christoffersen_test(realized_losses, var_forecasts, alpha)
  
  obs_uc  <- obs_kupiec$LRuc
  obs_ind <- obs_christ$LRind
  obs_cc  <- obs_christ$LRcc
  
  count_uc  <- 0L
  count_ind <- 0L
  count_cc  <- 0L
  
  for (i in seq_len(n_sim)) {
    sim_seq <- as.integer(runif(t_total) < p)
    
    x_sim <- sum(sim_seq)
    sim_uc <- lruc_stat(x_sim, t_total, p)
    
    v_today    <- sim_seq[-t_total]
    v_tomorrow <- sim_seq[-1]
    
    n00 <- sum(v_today == 0L & v_tomorrow == 0L)
    n01 <- sum(v_today == 0L & v_tomorrow == 1L)
    n10 <- sum(v_today == 1L & v_tomorrow == 0L)
    n11 <- sum(v_today == 1L & v_tomorrow == 1L)
    
    if ((n00 + n01) == 0L || (n10 + n11) == 0L) {
      sim_ind <- 0
      sim_cc  <- sim_uc
    } else {
      pi01 <- n01 / (n00 + n01)
      pi11 <- n11 / (n10 + n11)
      pi2  <- (n01 + n11) / (t_total - 1L)
      
      log_L0 <- xlogp(n00 + n10, 1 - pi2) + xlogp(n01 + n11, pi2)
      log_LA <- xlogp(n00, 1 - pi01) + xlogp(n01, pi01) +
        xlogp(n10, 1 - pi11) + xlogp(n11, pi11)
      
      sim_ind <- max(0, -2 * (log_L0 - log_LA))
      sim_cc  <- sim_uc + sim_ind
    }
    
    if (!is.na(sim_uc) && !is.na(obs_uc) && sim_uc >= obs_uc) count_uc <- count_uc + 1L
    if (!is.na(sim_ind) && !is.na(obs_ind) && sim_ind >= obs_ind) count_ind <- count_ind + 1L
    if (!is.na(sim_cc) && !is.na(obs_cc) && sim_cc >= obs_cc) count_cc <- count_cc + 1L
  }
  
  list(
    MC_p_uc  = count_uc / n_sim,
    MC_p_ind = if (is.na(obs_ind)) NA_real_ else count_ind / n_sim,
    MC_p_cc  = if (is.na(obs_cc)) NA_real_ else count_cc / n_sim
  )
}

# Acerbi-Szekely (2014) Z2: Z2 > -0.7 = green, -1.8 < Z2 <= -0.7 = yellow,
# Z2 <= -1.8 = red. Thresholds match the ones cited in Section 4.5.4.
acerbi_szekely_z2 <- function(realized_losses, var_forecasts, es_forecasts,
                              alpha = cfg$alpha_es) {
  ok <- !is.na(realized_losses) & !is.na(var_forecasts) & !is.na(es_forecasts)
  L   <- realized_losses[ok]
  VaR <- var_forecasts[ok]
  ES  <- es_forecasts[ok]
  T   <- length(L)
  
  if (T == 0L) return(list(Z2 = NA_real_, Classification = NA_character_))
  
  I_t <- as.numeric(L > VaR)
  Z2  <- 1 - sum((L * I_t) / ((1 - alpha) * ES)) / T
  
  zone <- if (Z2 <= -1.80)      "Red"
  else if (Z2 <= -0.70) "Yellow"
  else                  "Green"
  
  list(Z2 = Z2, Classification = zone)
}

# Basel 1996 thresholds rescaled from the canonical 250-day case via the
# Binomial CDF. Yellow entry = qbinom(0.95, T, p); Red entry = qbinom(0.9999).
basel_traffic_light <- function(t_total, violations, alpha = cfg$alpha_var) {
  if (t_total == 0L) return(NA_character_)
  p <- 1 - alpha
  yellow_entry <- qbinom(0.95,   t_total, p)
  red_entry    <- qbinom(0.9999, t_total, p)
  
  if (violations >= red_entry)         "Red"
  else if (violations >= yellow_entry) "Yellow"
  else                                 "Green"
}