# --- STEP 3: STRUCTURAL BREAKS AND MARKET ASYMMETRY ---
# This script analyzes the temporal stability of beta and its behavior in different market conditions.

library(quantmod)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(strucchange)
library(sandwich)
library(lmtest)

# Load data processed in previous steps
load("capm_data.RData")
# PERIOD DEFINITION ------------------------------------------------------------

# The sample is split into two parts:
# Period 1: Feb 2019 - Dec 2023 (59 months)
# Period 2: Jan 2024 - Dec 2025 (24 months)

breakpoint <- 59
n_total    <- length(mkt_prem)

cat("Total market premium length:", n_total, "months\n")
cat("Period 1: Indices 1 to", breakpoint, "\n")
cat("Period 2: Indices", breakpoint + 1, "to", n_total, "\n")


# 1. CHOW TESTS FOR STRUCTURAL STABILITY ---------------------------------------
chow_results <- data.frame()

for (tk in names(returns_list)) {
  
  stock_prem <- returns_list[[tk]] - rf_ret
  n          <- min(length(stock_prem), n_total)
  y          <- stock_prem[1:n]
  x          <- mkt_prem[1:n]
  valid      <- complete.cases(x, y)
  
  n1 <- sum(valid[1:breakpoint])
  n2 <- sum(valid[(breakpoint + 1):n])
  if (n1 < 20 || n2 < 20) next # Minimum observation requirement
  
  y_v      <- y[valid]
  x_v      <- x[valid]
  bp_valid <- sum(valid[1:breakpoint])
  
  model_full <- lm(y_v ~ x_v)
  model_1    <- lm(y_v[1:bp_valid] ~ x_v[1:bp_valid])
  model_2    <- lm(y_v[(bp_valid+1):length(y_v)] ~ x_v[(bp_valid+1):length(x_v)])
  
  RSS_full <- sum(residuals(model_full)^2)
  RSS_1    <- sum(residuals(model_1)^2)
  RSS_2    <- sum(residuals(model_2)^2)
  k        <- 2 
  n_obs    <- length(y_v)
  
  # F-statistic: Comparing RSS of full model vs. pooled sub-models
  F_chow <- ((RSS_full - (RSS_1 + RSS_2)) / k) /
    ((RSS_1 + RSS_2) / (n_obs - 2 * k))
  p_chow <- pf(F_chow, df1 = k, df2 = n_obs - 2 * k, lower.tail = FALSE)
  
  # Forecast Test: Applying Period 1 estimates to Period 2
  df_v   <- data.frame(y = y_v, x = x_v)
  sc     <- sctest(y ~ x, data = df_v, type = "Chow", point = bp_valid)
  F_pred <- sc$statistic
  p_pred <- sc$p.value
  
  beta_1 <- coef(model_1)[2]
  beta_2 <- coef(model_2)[2]
  
  chow_results <- rbind(chow_results, data.frame(
    ticker     = tk,
    beta_1     = round(beta_1, 4),
    beta_2     = round(beta_2, 4),
    delta_beta = round(beta_2 - beta_1, 4),
    F_chow     = round(F_chow, 3),
    p_chow     = round(p_chow, 4),
    F_pred     = round(F_pred, 3),
    p_pred     = round(p_pred, 4),
    chow_sig   = ifelse(p_chow < 0.05, "Break (5%)", "No Break"),
    pred_sig   = ifelse(p_pred < 0.05, "Break (5%)", "No Break")
  ))
}

row.names(chow_results) <- NULL

# Overview of observation counts
chow_n <- data.frame()

for (tk in names(returns_list)) {
  stock_prem <- returns_list[[tk]] - rf_ret
  n          <- min(length(stock_prem), length(mkt_prem))
  y          <- stock_prem[1:n]
  x          <- mkt_prem[1:n]
  valid      <- complete.cases(x, y)
  
  n1 <- sum(valid[1:breakpoint])
  n2 <- sum(valid[(breakpoint + 1):n])
  
  chow_n <- rbind(chow_n, data.frame(
    ticker   = tk,
    n1       = n1,
    n2       = n2,
    n_total  = n1 + n2,
    status   = ifelse(n1 < 20 | n2 < 20, "Skipped", "Included")
  ))
}

cat("\n--- Observation Counts per Period ---\n")
print(chow_n %>% arrange(n1) %>% head(5), row.names = FALSE)

cat("\n=== Chow Test Results ===\n")
print(chow_results %>% arrange(p_pred) %>%
        select(ticker, beta_1, beta_2, delta_beta,
               p_chow, chow_sig, p_pred, pred_sig))


# Plotting Beta change sorted by magnitude
p_chow_delta <- chow_results %>%
  arrange(delta_beta) %>%
  mutate(ticker = factor(ticker, levels = ticker)) %>%
  ggplot(aes(x = ticker, y = delta_beta, fill = chow_sig)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
  scale_fill_manual(values = c("Break (5%)" = "#CC0000",
                                "No Break"  = "#003087")) +
  labs(title    = "Change in Beta Coefficient (β₂ − β₁)",
       subtitle = "Red = Statistically significant structural break (5%)",
       x = NULL, y = "Δβ = β(2024–2025) − β(2019–2023)",
       fill = "Status") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

print(p_chow_delta)
ggsave("chow_delta_beta.png", p_chow_delta, width = 12, height = 5, dpi = 300)

# 2. SML ANALYSIS ACROSS TWO PERIODS -------------------------------------------

# Define indices
idx_1 <- 1:breakpoint
idx_2 <- (breakpoint + 1):n_total

# Statistics for Period 1 and 2
stats_1 <- data.frame()
stats_2 <- data.frame()

for (p in 1:10) {
  tickers_p <- capm_results$ticker[capm_results$portfolio == p]
  
  # A) Average portfolio excess return in Period 1
  port_rets_list <- lapply(tickers_p, function(tk) returns_list[[tk]][idx_1])
  ret_matrix_1   <- do.call(cbind, port_rets_list)
  avg_ret_p1     <- rowMeans(ret_matrix_1, na.rm = TRUE)
  avg_prem_p1    <- mean(avg_ret_p1 - rf_ret[idx_1], na.rm = TRUE)
  
  # B) Average portfolio Beta in Period 1
  beta_vals_p1 <- sapply(tickers_p, function(tk) {
    y <- (returns_list[[tk]] - rf_ret)[idx_1]
    x <- mkt_prem[idx_1]
    v <- complete.cases(y, x)
    if (sum(v) < 20) return(NA)
    coef(lm(y[v] ~ x[v]))[2]
  })
  
  stats_1 <- rbind(stats_1, data.frame(
    portfolio = p, 
    beta_p    = mean(beta_vals_p1, na.rm = TRUE), 
    avg_prem  = avg_prem_p1
  ))
}

for (p in 1:10) {
  tickers_p <- capm_results$ticker[capm_results$portfolio == p]
  
  # A) Average portfolio excess return in Period 2
  port_rets_list <- lapply(tickers_p, function(tk) returns_list[[tk]][idx_2])
  ret_matrix_2   <- do.call(cbind, port_rets_list)
  avg_ret_p2     <- rowMeans(ret_matrix_2, na.rm = TRUE)
  avg_prem_p2    <- mean(avg_ret_p2 - rf_ret[idx_2], na.rm = TRUE)
  
  # B) Average portfolio Beta in Period 2
  beta_vals_p2 <- sapply(tickers_p, function(tk) {
    y <- (returns_list[[tk]] - rf_ret)[idx_2]
    x <- mkt_prem[idx_2]
    v <- complete.cases(y, x)
    if (sum(v) < 20) return(NA)
    coef(lm(y[v] ~ x[v]))[2]
  })
  
  stats_2 <- rbind(stats_2, data.frame(
    portfolio = p, 
    beta_p    = mean(beta_vals_p2, na.rm = TRUE), 
    avg_prem  = avg_prem_p2
  ))
}

# SML Models for both periods
sml_model_1 <- lm(avg_prem ~ beta_p, data = stats_1)
sml_rob_1   <- coeftest(sml_model_1, vcov = vcovHAC(sml_model_1))
print(sml_rob_1)

sml_model_2 <- lm(avg_prem ~ beta_p, data = stats_2)
sml_rob_2   <- coeftest(sml_model_2, vcov = vcovHAC(sml_model_2))
print(sml_rob_2)

# Theoretical market premiums
avg_mkt_1 <- mean(mkt_prem[idx_1], na.rm = TRUE)
avg_mkt_2 <- mean(mkt_prem[idx_2], na.rm = TRUE)

x_lims <- range(c(stats_1$beta_p, stats_2$beta_p), na.rm = TRUE)
y_lims <- range(c(stats_1$avg_prem, stats_2$avg_prem), na.rm = TRUE)

x_expand <- c(x_lims[1] * 0.9, x_lims[2] * 1.1)
y_expand <- c(min(0, y_lims[1] * 1.1), y_lims[2] * 1.1)

# Period 1 Plot
plot_1 <- ggplot(stats_1, aes(x = beta_p, y = avg_prem)) +
  geom_point(size = 3, color = "#003087") +
  geom_text(aes(label = paste0("P", portfolio)), vjust = -0.8, size = 3, color = "#003087") +
  geom_smooth(method = "lm", se = TRUE, color = "#CC0000", linewidth = 1) +
  geom_abline(intercept = 0, slope = avg_mkt_1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  scale_x_continuous(limits = x_expand) +
  scale_y_continuous(limits = y_expand) +
  labs(title = "SML: Period 1 (Stable)", x = "Beta", y = "Return")

# Period 2 Plot
plot_2 <- ggplot(stats_2, aes(x = beta_p, y = avg_prem)) +
  geom_point(size = 3, color = "#003087") +
  geom_text(aes(label = paste0("P", portfolio)), vjust = -0.8, size = 3, color = "#003087") +
  geom_smooth(method = "lm", se = TRUE, color = "#CC0000", linewidth = 1) +
  geom_abline(intercept = 0, slope = avg_mkt_2, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  scale_x_continuous(limits = x_expand) +
  scale_y_continuous(limits = y_expand) +
  labs(title = "SML: Period 2 (Last 24m)", x = "Beta", y = "Return")

final_plot_bp <- grid.arrange(plot_1, plot_2, ncol = 2)
ggsave("SML_comparison.png", final_plot_bp, width = 12, height = 6, dpi = 300)


# 3. ASYMMETRIC BETA ANALYSIS --------------------------------------------------
D_up   <- as.numeric(mkt_prem > 0)
D_down <- 1 - D_up

cat(sprintf("\nUp-Market Months (Growth):   %d (%.1f%%)\n", sum(D_up),   100 * mean(D_up)))
cat(sprintf("Down-Market Months (Decline): %d (%.1f%%)\n",   sum(D_down), 100 * mean(D_down)))

asym_results <- data.frame()

for (tk in names(returns_list)) {
  stock_prem <- returns_list[[tk]] - rf_ret
  n     <- min(length(stock_prem), n_total)
  y     <- stock_prem[1:n]
  x     <- mkt_prem[1:n]
  d_up  <- D_up[1:n]
  d_dn  <- D_down[1:n]
  valid <- complete.cases(y, x)
  if (sum(valid) < 30) next
  
  x_up <- x * d_up
  x_dn <- x * d_dn
  
  model_asym <- lm(y[valid] ~ x_up[valid] + x_dn[valid])
  cf         <- coef(model_asym)
  beta_up    <- cf[2]
  beta_dn    <- cf[3]
  
  model_sym <- lm(y[valid] ~ x[valid])
  f_test    <- anova(model_sym, model_asym)
  p_asym    <- f_test$`Pr(>F)`[2]
  
  asym_results <- rbind(asym_results, data.frame(
    ticker     = tk,
    beta_up    = round(beta_up, 4),
    beta_dn    = round(beta_dn, 4),
    delta_beta = round(beta_up - beta_dn, 4),
    p_asym     = round(p_asym, 4),
    asym_sig   = ifelse(p_asym < 0.05, "Asymmetry (5%)", "Symmetry")
  ))
}

row.names(asym_results) <- NULL

# Average beta values for comparison
colMeans(asym_results[, c("beta_up", "beta_dn", "delta_beta")], na.rm = TRUE)

cat(sprintf("\nNumber of stocks with asymmetric beta: %d out of %d (%.1f%%)\n",
            sum(asym_results$p_asym < 0.05),
            nrow(asym_results),
            100 * mean(asym_results$p_asym < 0.05)))

cat("\n=== Asymmetry Test Results ===\n")
print(asym_results %>% arrange(p_asym) %>%
        select(ticker, beta_up, beta_dn, delta_beta, p_asym, asym_sig))

# Asymmetric Beta Plot
p_asym_scatter <- ggplot(asym_results, aes(x = beta_dn, y = beta_up)) +
  geom_point(aes(color = asym_sig), size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Asymmetry (5%)" = "#CC0000",
                                "Symmetry"       = "#003087")) +
  labs(title    = "Asymmetric Beta – Up vs. Down Market",
       subtitle = "Above Diagonal: β⁺ > β⁻ (Asset is more sensitive to growth)",
       x = "Beta in Down-Market (β⁻)", y = "Beta in Up-Market (β⁺)",
       color = "Status") +
  theme_minimal(base_size = 11)

print(p_asym_scatter)
ggsave("asymmetric_beta.png", p_asym_scatter, width = 8, height = 6, dpi = 300)

# 4. SML FOR UP AND DOWN MARKET ------------------------------------------------

# Index preparation
n_mkt <- length(mkt_prem)
idx_up <- which(D_up == 1)
idx_dn <- which(D_down == 1)

# UP-MARKET SML
stats_up <- data.frame()

for (p in 1:10) {
  tickers_p <- capm_results$ticker[capm_results$portfolio == p]
  
  ret_matrix <- do.call(cbind, lapply(tickers_p, function(tk) {
    r <- returns_list[[tk]]
    out <- rep(NA, n_mkt)
    n <- min(length(r), n_mkt)
    out[1:n] <- r[1:n]
    out[idx_up]
  }))
  
  port_ret  <- if (is.null(dim(ret_matrix))) ret_matrix else rowMeans(ret_matrix, na.rm = TRUE)
  port_prem <- port_ret - rf_ret[idx_up]
  valid     <- complete.cases(port_prem)
  
  beta_vals_raw <- sapply(tickers_p, function(tk) {
    val <- asym_results$beta_up[asym_results$ticker == tk]
    if(length(val) == 0) return(NA) else return(as.numeric(val))
  })
  beta_clean <- as.numeric(unname(unlist(beta_vals_raw)))
  
  stats_up <- rbind(stats_up, data.frame(
    portfolio = p,
    beta_p    = mean(beta_clean, na.rm = TRUE),
    avg_prem  = mean(port_prem[valid], na.rm = TRUE)
  ))
}
# SML Model for Up-Market
sml_m_up  <- lm(avg_prem ~ beta_p, data = stats_up)
sml_rob_up <- coeftest(sml_m_up, vcov = vcovHAC(sml_m_up))
avg_mkt_up <- mean(mkt_prem[idx_up], na.rm = TRUE)

gamma0_up <- sml_rob_up[1, 1]
gamma1_up <- sml_rob_up[2, 1]
p_g0_up   <- sml_rob_up[1, 4]
p_g1_up   <- sml_rob_up[2, 4]

# DOWN-MARKET SML
stats_dn <- data.frame()

for (p in 1:10) {
  tickers_p <- capm_results$ticker[capm_results$portfolio == p]
  
  ret_matrix <- do.call(cbind, lapply(tickers_p, function(tk) {
    r <- returns_list[[tk]]
    out <- rep(NA, n_mkt)
    n <- min(length(r), n_mkt)
    out[1:n] <- r[1:n]
    out[idx_dn]
  }))
  
  port_ret  <- if (is.null(dim(ret_matrix))) ret_matrix else rowMeans(ret_matrix, na.rm = TRUE)
  port_prem <- port_ret - rf_ret[idx_dn]
  valid     <- complete.cases(port_prem)
  
  beta_vals_raw <- sapply(tickers_p, function(tk) {
    val <- asym_results$beta_dn[asym_results$ticker == tk]
    if(length(val) == 0) return(NA) else return(as.numeric(val))
  })
  beta_clean <- as.numeric(unname(unlist(beta_vals_raw)))
  
  stats_dn <- rbind(stats_dn, data.frame(
    portfolio = p,
    beta_p    = mean(beta_clean, na.rm = TRUE),
    avg_prem  = mean(port_prem[valid], na.rm = TRUE)
  ))
}

# SML Model for Down-Market
sml_m_dn  <- lm(avg_prem ~ beta_p, data = stats_dn)
sml_rob_dn <- coeftest(sml_m_dn, vcov = vcovHAC(sml_m_dn))
avg_mkt_dn <- mean(mkt_prem[idx_dn], na.rm = TRUE)

gamma0_dn <- sml_rob_dn[1, 1]
gamma1_dn <- sml_rob_dn[2, 1]
p_g0_dn   <- sml_rob_dn[1, 4]
p_g1_dn   <- sml_rob_dn[2, 4]

cat("\n--- SML: Up-Market Period (Robust Newey-West) ---\n")
cat(sprintf("γ₀ = %.5f (p = %.4f)\n", gamma0_up, p_g0_up))
cat(sprintf("γ₁ = %.5f (p = %.5f)\n", gamma1_up, p_g1_up))
cat(sprintf("Market Risk Premium = %.5f\n", avg_mkt_up))

cat("\n--- SML: Down-Market Period (Robust Newey-West) ---\n")
cat(sprintf("γ₀ = %.5f (p = %.4f)\n", gamma0_dn, p_g0_dn))
cat(sprintf("γ₁ = %.5f (p = %.4f)\n", gamma1_dn, p_g1_dn))
cat(sprintf("Market Risk Premium = %.5f\n", avg_mkt_dn))
# Axis range preparation
x_min <- min(c(stats_up$beta_p, stats_dn$beta_p), na.rm = TRUE) * 0.9
x_max <- max(c(stats_up$beta_p, stats_dn$beta_p), na.rm = TRUE) * 1.1
y_abs_max <- max(abs(c(stats_up$avg_prem, stats_dn$avg_prem)), na.rm = TRUE) * 1.2

common_x <- scale_x_continuous(limits = c(x_min, x_max))
common_y <- scale_y_continuous(limits = c(-y_abs_max, y_abs_max))

# Mirror Chart: Up Market
plot_up <- ggplot(stats_up, aes(x = beta_p, y = avg_prem)) +
  geom_point(size = 3, color = "#228B22") + 
  geom_text(aes(label = paste0("P", portfolio)), vjust = -1, size = 3.5, color = "#228B22") +
  geom_smooth(method = "lm", formula = y ~ x, color = "darkgreen", se = TRUE) +
  geom_abline(intercept = 0, slope = avg_mkt_up, linetype = "dashed", color = "grey60") +
  common_x + common_y +
  labs(title    = "SML: Up-Market (Robust)",
       subtitle = sprintf("γ₀=%.4f (p=%.3f), γ₁=%.4f (p=%.3f)", 
                          gamma0_up, p_g0_up, gamma1_up, p_g1_up),
       x = "Beta (Up)", y = "Mean Excess Return") +
  theme_minimal()

# Mirror Chart: Down Market
plot_dn <- ggplot(stats_dn, aes(x = beta_p, y = avg_prem)) +
  geom_point(size = 3, color = "#CC0000") + 
  geom_text(aes(label = paste0("P", portfolio)), vjust = 1.5, size = 3.5, color = "#CC0000") +
  geom_smooth(method = "lm", formula = y ~ x, color = "darkred", se = TRUE) +
  geom_abline(intercept = 0, slope = avg_mkt_dn, linetype = "dashed", color = "grey60") +
  common_x + common_y +
  labs(title    = "SML: Down-Market (Robust)",
       subtitle = sprintf("γ₀=%.4f (p=%.3f), γ₁=%.4f (p=%.3f)", 
                          gamma0_dn, p_g0_dn, gamma1_dn, p_g1_dn),
       x = "Beta (Down)", y = "Mean Excess Return") +
  theme_minimal()

final_plot_mirror <- grid.arrange(plot_up, plot_dn, ncol = 2)
ggsave("SML_up_vs_down.png", final_plot_mirror, width = 12, height = 6, dpi = 300)

# 5. DATA EXPORT ---------------------------------------------------------------
save(capm_results, returns_list, mkt_prem, rf_ret, 
     dates_ret, n100_vec, dates_index, rf_vec, 
     portfolio_stats, chow_results, asym_results, stats_up, stats_dn,
     file = "capm_data.RData")
