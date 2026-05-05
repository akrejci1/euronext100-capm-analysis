# --- STEP 1: DATA ACQUISITION & BETA ESTIMATION ---
# Part of the Euronext 100 CAPM Analysis Project

library(quantmod)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(strucchange)
library(sandwich)
library(lmtest)

# 1. TICKERS SELECTION ---------------------------------------------------------
# 100 constituents according to the official Euronext list
tickers_n100 <- c(
  # --- Netherlands / Amsterdam (21) ---
  "ABN.AS", "ADYEN.AS", "AGN.AS", "AD.AS", "AKZA.AS", "MT.AS", "ASM.AS", 
  "ASML.AS", "ASRNL.AS", "DSFIR.AS", "EXO.AS", "HEIA.AS", "INGA.AS", 
  "JDEP.AS", "KPN.AS", "NN.AS", "PHIA.AS", "PRX.AS", "SHELL.AS", "UMG.AS", "WKL.AS",
  
  # --- Belgium / Brussels (5) ---
  "ABI.BR", "AGS.BR", "ARGX.BR", "KBC.BR", "UCB.BR",
  
  # --- Ireland / Dublin (5) ---
  "A5G.IR", "BIRG.IR", "KRZ.IR", "KRX.IR", "RYA.IR",
  
  # --- Portugal / Lisbon (3) ---
  "EDP.LS", "GALP.LS", "JMT.LS",
  
  # --- Italy / Milan (22) ---
  "BMED.MI", "BAMI.MI", "BPE.MI", "ENEL.MI", "ENI.MI", "RACE.MI", "FBK.MI", 
  "G.MI", "ISP.MI", "INW.MI", "LDO.MI", "MB.MI", "MONC.MI", "PST.MI", "PRY.MI", 
  "REC.MI", "SRG.MI", "STLAM.MI", "TEN.MI", "TRN.MI", "UCG.MI", "UNI.MI",

  # --- France / Paris (38) ---
  "AC.PA", "ADP.PA", "AI.PA", "AIR.PA", "CS.PA", "BIM.PA", "BNP.PA", "EN.PA", 
  "BVI.PA", "CAP.PA", "ACA.PA", "BN.PA", "DSY.PA", "FGR.PA", "ENGI.PA", "EL.PA", 
  "ERF.PA", "ENX.PA", "IPN.PA", "KER.PA", "LR.PA", "MC.PA", "ML.PA", "ORA.PA", 
  "RI.PA", "PUB.PA", "RNO.PA", "SAF.PA", "SGO.PA", "SAN.PA", "SU.PA", "GLE.PA", 
  "STMPA.PA", "HO.PA", "TTE.PA", "URW.PA", "VIE.PA", "DG.PA",
  
  # --- Norway / Oslo (6) ---
  "AKRBP.OL", "DNB.OL", "EQNR.OL", "KOG.OL", "NHY.OL", "TEL.OL"
)

cat("Total tickers selected:", length(tickers_n100), "\n")

# 2. DATA DOWNLOAD -------------------------------------------------------------
start_date <- "2019-01-01"
end_date   <- "2025-12-31"

# Market Index: Euronext 100 (^N100) from Yahoo Finance
cat("\nDownloading Euronext 100 index (^N100)...\n")
n100_raw     <- getSymbols("^N100", src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
n100_monthly <- to.monthly(Ad(n100_raw), indexAt = "lastof", OHLC = FALSE)
n100_monthly <- na.omit(n100_monthly, na.rm = FALSE)

# Clean NAs after interpolation
first_valid  <- min(which(!is.na(as.numeric(n100_monthly))))
last_valid   <- max(which(!is.na(as.numeric(n100_monthly))))
n100_monthly <- n100_monthly[first_valid:last_valid]

dates_index  <- index(n100_monthly)
n100_vec     <- as.numeric(n100_monthly)
cat("  Successfully loaded:", length(n100_vec), "months\n")

# Risk-Free Rate: ECB Deposit Facility Rate (FRED: ECBDFR)[cite: 5]
cat("Downloading Risk-Free Rate (FRED: ECBDFR)...\n")
rf_vec <- tryCatch({
  rf_raw    <- getSymbols("ECBDFR", src = "FRED", from = start_date, to = end_date, auto.assign = FALSE)
  rf_m      <- to.monthly(rf_raw, indexAt = "lastof", OHLC = FALSE)
  rf_approx <- approx(as.numeric(index(rf_m)), as.numeric(rf_m), xout = as.numeric(dates_index), rule = 2)$y
  cat("  Successfully loaded\n")
  rf_approx / 100 / 12
}, error = function(e) {
  cat("  FRED unavailable – setting rf = 0\n")
  rep(0, length(dates_index))
})

# Stock Prices from Yahoo Finance[cite: 5]
cat("\nDownloading stock prices...\n")
prices_list    <- list()
failed_tickers <- character(0)

for (tk in tickers_n100) {
  tryCatch({
    dat     <- getSymbols(tk, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
    monthly <- to.monthly(Ad(dat), indexAt = "lastof", OHLC = FALSE)
    aligned <- merge(n100_monthly, monthly, join = "left")[, 2]
    prices_list[[tk]] <- as.numeric(aligned)
    cat("  Loaded:", tk, "\n")
  }, error = function(e) {
    cat("  FAILED:", tk, "\n")
    failed_tickers <<- c(failed_tickers, tk)
  })
}

cat("\n--- Download Summary ---\n")
cat("Successfully downloaded:", length(prices_list), "tickers\n")
if (length(failed_tickers) > 0) cat("Failed tickers:", paste(failed_tickers, collapse = ", "), "\n")

# 3. CALCULATE MONTHLY LOGARITHMIC RETURNS -------------------------------------
log_ret <- function(p) {
  r <- diff(log(p))
  r[is.infinite(r) | is.nan(r)] <- NA
  r
}

T            <- length(n100_vec)
rm_vec       <- log_ret(n100_vec)
rf_ret       <- tail(rf_vec, T - 1)
mkt_prem     <- rm_vec - rf_ret
dates_ret    <- tail(dates_index, T - 1)
returns_list <- lapply(prices_list, log_ret)

# Data availability check[cite: 5]
data_overview <- do.call(rbind, lapply(names(returns_list), function(tk) {
  r       <- returns_list[[tk]]
  valid_i <- which(!is.na(r))
  data.frame(
    ticker    = tk,
    n_valid   = length(valid_i),
    n_missing = sum(is.na(r)),
    pct_valid = round(100 * length(valid_i) / length(r), 1),
    first_obs = as.character(dates_ret[min(valid_i)]),
    last_obs  = as.character(dates_ret[max(valid_i)])
  )
}))

data_overview <- data_overview %>% arrange(n_valid)
cat("\n--- Data Availability (sorted by number of valid months) ---\n")
print(data_overview, row.names = FALSE)

# 4. VISUALIZATION OF INPUTS ---------------------------------------------------
# N100 Index Plot
N100g <- ggplot(data.frame(date = dates_index, N100 = n100_vec), aes(x = date, y = N100)) +
  geom_line(color = "#003087", linewidth = 1) +
  labs(title = "Euronext 100 Index Performance (2019–2025)", x = NULL, y = "Index Value") +
  theme_minimal(base_size = 11)
print(N100g)
ggsave("N100_performance.png", N100g, width = 9, height = 4, dpi = 300)

# ECB Rate Plot
ECB_rate <- ggplot(data.frame(date = dates_index, rf = rf_vec * 100 * 12), aes(x = date, y = rf)) +
  geom_line(color = "#CC0000", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "ECB Deposit Facility Rate (% p.a.)", x = NULL, y = "Rate (%)") +
  theme_minimal(base_size = 11)
print(ECB_rate)
ggsave("ECB_rate_performance.png", ECB_rate, width = 9, height = 4, dpi = 150)

# Mean Returns Distribution
mean_rets <- sapply(returns_list, mean, na.rm = TRUE)
dist_returns <- ggplot(data.frame(mean_ret = mean_rets), aes(x = mean_ret)) +
  geom_histogram(bins = 20, fill = "#003087", color = "white", alpha = 0.85) +
  geom_vline(xintercept = mean(mkt_prem, na.rm = TRUE), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of Mean Monthly Returns – Euronext 100",
       subtitle = "Red line = Mean market risk premium (^N100), 2019–2025",
       x = "Mean Monthly Return", y = "Frequency") +
  theme_minimal(base_size = 11)
print(dist_returns)
ggsave("returns_distribution.png", dist_returns, width = 8, height = 4, dpi = 300)

# 5. BETA COEFFICIENT ESTIMATION -----------------------------------------------
capm_results <- data.frame()
char_data    <- list()

for (tk in names(returns_list)) {
  stock_prem <- returns_list[[tk]] - rf_ret
  n     <- min(length(stock_prem), length(mkt_prem))
  y     <- stock_prem[1:n]
  x     <- mkt_prem[1:n]
  valid <- complete.cases(x, y)
  if (sum(valid) < 30) { cat("Insufficient data:", tk, "\n"); next }
  
  model <- lm(y[valid] ~ x[valid])
  s_robust <- coeftest(model, vcov = vcovHAC(model))
  cf <- s_robust
  
  capm_results <- rbind(capm_results, data.frame(
    ticker     = tk,
    alpha      = cf[1, 1],
    beta       = cf[2, 1],
    se_alpha   = cf[1, 2],
    se_beta    = cf[2, 2],
    t_alpha    = cf[1, 3],
    p_alpha    = cf[1, 4],
    r_squared  = summary(model)$r.squared,
    sigma2_eps = var(residuals(model))
  ))
  char_data[[tk]] <- list(x = x[valid], y = y[valid])
}

capm_results <- capm_results %>%
  arrange(beta) %>%
  mutate(
    alpha_sig = ifelse(p_alpha < 0.05, "Reject H0 (α≠0)", "Fail to reject H0 (α=0)"),
    beta_type = ifelse(beta >= 1, "Aggressive (β≥1)", "Defensive (β<1)")[cite: 5]
  )

print(capm_results %>% arrange(p_alpha) %>% 
        select(ticker, alpha, beta, r_squared, p_alpha, alpha_sig), digits = 4)

write.csv(capm_results, "capm_beta_results.csv", row.names = FALSE)

# 7. CHARACTERISTIC LINES VISUALIZATION ----------------------------------------
sel_idx   <- round(seq(1, nrow(capm_results), length.out = 9))
selected  <- capm_results$ticker[sel_idx]
plot_list <- list()

for (tk in selected) {
  row     <- capm_results[capm_results$ticker == tk, ]
  df_plot <- data.frame(x = char_data[[tk]]$x, y = char_data[[tk]]$y)
  plot_list[[tk]] <- ggplot(df_plot, aes(x = x, y = y)) +
    geom_point(alpha = 0.45, size = 1.5, color = "#003087") +
    geom_smooth(method = "lm", se = TRUE, color = "#CC0000", linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    labs(title    = tk,
         subtitle = sprintf("β=%.3f | α=%.4f | R²=%.3f", row$beta, row$alpha, row$r_squared),
         x = expression(r[m]-r[f]), y = expression(r[j]-r[f])) +
    theme_minimal(base_size = 8)
}

final_grid <- do.call(grid.arrange, c(plot_list, ncol = 3, 
          top = "CAPM Characteristic Lines – Selection of 9 Assets (Euronext 100)"))

ggsave("characteristic_lines.png", plot = final_grid, width = 12, height = 10, dpi = 300)

# 8. DATA EXPORT ---------------------------------------------------------------
save(capm_results, returns_list, mkt_prem, rf_ret,
     dates_ret, n100_vec, dates_index, rf_vec,
     file = "capm_data.RData")
