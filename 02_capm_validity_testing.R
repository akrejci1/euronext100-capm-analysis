# --- STEP 2: CAPM VALIDITY TESTING ---
# This script tests the Security Market Line (SML) and performs Fama-MacBeth regressions.

library(quantmod)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(strucchange)
library(sandwich)
library(lmtest)

# Load data processed in Step 1
load("capm_data.RData")

# 1. PORTFOLIO CONSTRUCTION (10 BETA-SORTED PORTFOLIOS) -----------------------

n_portfolios <- 10

# Rank assets by estimated beta and divide into deciles
capm_results <- capm_results %>%
  mutate(portfolio = ntile(beta, n_portfolios))

cat("Distribution of assets into portfolios:\n")
print(capm_results %>%
        group_by(portfolio) %>%
        summarise(n = n(), beta_min = min(beta), beta_max = max(beta),
                  beta_mean = mean(beta)) %>%
        as.data.frame(), digits = 3)

# 2. CALCULATE PORTFOLIO RETURNS AND BETAS ------------------------------------

portfolio_stats <- data.frame()

for (p in 1:n_portfolios) {
  tickers_p  <- capm_results$ticker[capm_results$portfolio == p]
  n_mkt      <- length(mkt_prem)
  
  # Align monthly returns for all assets in the portfolio
  ret_matrix <- do.call(cbind, lapply(tickers_p, function(tk) {
    r <- returns_list[[tk]]
    n <- min(length(r), n_mkt)
    out <- rep(NA, n_mkt)
    out[1:n] <- r[1:n]
    out
  }))
  
  # Calculate equal-weighted portfolio return
  port_ret <- rowMeans(ret_matrix, na.rm = TRUE)
  
  # Calculate portfolio excess return
  port_prem <- port_ret - rf_ret
  valid     <- complete.cases(port_prem, mkt_prem)
  
  # Calculate aggregate portfolio statistics
  beta_p    <- mean(capm_results$beta[capm_results$portfolio == p])
  sigma2_p  <- mean(capm_results$sigma2_eps[capm_results$portfolio == p])
  
  portfolio_stats <- rbind(portfolio_stats, data.frame(
    portfolio  = p,
    n_assets   = length(tickers_p),
    beta_p     = beta_p,
    avg_prem   = mean(port_prem[valid], na.rm = TRUE),
    sigma2_eps = sigma2_p,
    tickers    = paste(tickers_p, collapse = ", ")
  ))
}

cat("\nPortfolio Results Summary:\n")
print(portfolio_stats %>% select(portfolio, n_assets, beta_p, avg_prem, sigma2_eps),
      digits = 4)

# 3. SECURITY MARKET LINE (SML) ESTIMATION ------------------------------------

avg_mkt_prem <- mean(mkt_prem, na.rm = TRUE)
sml_model    <- lm(avg_prem ~ beta_p, data = portfolio_stats)

# Calculate robust standard errors (HAC)
sml_rob    <- coeftest(sml_model, vcov = vcovHAC(sml_model))
gamma0_sml <- sml_rob[1, 1]
gamma1_sml <- sml_rob[2, 1]

cat("\n=== Security Market Line – Cross-sectional Regression (Robust SE) ===\n")
print(sml_rob)
cat(sprintf("\nEstimate γ₀ = %.5f (p = %.4f)\n", gamma0_sml, sml_rob[1, 4]))
cat(sprintf("Estimate γ₁ = %.5f (p = %.4f)\n",   gamma1_sml, sml_rob[2, 4]))
cat(sprintf("Observed Market Risk Premium (rm-rf) = %.5f\n", avg_mkt_prem))

# Interpretation of the intercept (γ₀)
if (sml_rob[1, 4] > 0.05) {
  cat("  γ₀: Fail to reject H0 (γ₀=0) → Consistent with CAPM\n")
} else {
  cat("  γ₀: Reject H0 (γ₀=0) → Inconsistent with CAPM\n")
}

# Interpretation of the slope (γ₁)
if (abs(gamma1_sml - avg_mkt_prem) < 2 * sml_rob[2, 2]) {
  cat("  γ₁: Close to observed market premium → Consistent with CAPM\n")
} else {
  cat("  γ₁: Significant deviation from market premium → Inconsistent with CAPM\n")
}

# SML Visualization
p_sml <- ggplot(portfolio_stats, aes(x = beta_p, y = avg_prem)) +
  geom_point(size = 3, color = "#003087") +
  geom_text(aes(label = paste0("P", portfolio)),
            vjust = -0.8, size = 3, color = "#003087") +
  geom_smooth(method = "lm", se = TRUE, color = "#CC0000", linewidth = 1) +
  geom_abline(intercept = 0, slope = avg_mkt_prem,
              linetype = "dashed", color = "grey50", linewidth = 0.8) +
  labs(title = "Security Market Line – Euronext 100",
       subtitle = sprintf("γ₀=%.4f (p=%.3f), γ₁=%.4f (p=%.3f) | Grey = Theoretical SML",
                          gamma0_sml, sml_rob[1, 4], gamma1_sml, sml_rob[2, 4]),
       x = "Portfolio Beta (β)",
       y = "Mean Monthly Risk Premium (r_p - r_f)") +
  theme_minimal(base_size = 11)

print(p_sml)
ggsave("security_market_line.png", p_sml, width = 9, height = 5, dpi = 300)

# 4. FAMA-MACBETH REGRESSION (1973) -------------------------------------------

n_mkt            <- length(mkt_prem)
beta_vec         <- portfolio_stats$beta_p
beta2_vec        <- beta_vec^2
sigma2_vec       <- portfolio_stats$sigma2_eps

# Build portfolio monthly excess premium matrix (T × 10)
port_prem_matrix <- matrix(NA, nrow = n_mkt, ncol = n_portfolios)
for (p in 1:n_portfolios) {
  tickers_p  <- capm_results$ticker[capm_results$portfolio == p]
  ret_matrix <- do.call(cbind, lapply(tickers_p, function(tk) {
    r   <- returns_list[[tk]]
    out <- rep(NA, n_mkt)
    n   <- min(length(r), n_mkt)
    out[1:n] <- r[1:n]
    out
  }))
  port_ret <- if (is.null(dim(ret_matrix))) ret_matrix else rowMeans(ret_matrix, na.rm = TRUE)
  port_prem_matrix[, p] <- port_ret - rf_ret
}

# Perform monthly cross-sectional regressions
gammas_fm <- data.frame(
  t      = 1:n_mkt,
  gamma0 = NA_real_,
  gamma1 = NA_real_,
  gamma2 = NA_real_,
  gamma3 = NA_real_
)

for (t in 1:n_mkt) {
  y_t   <- port_prem_matrix[t, ]
  valid <- complete.cases(y_t, beta_vec, beta2_vec, sigma2_vec)
  if (sum(valid) < 6) next
  
  # Estimate using White HC3 robust errors
  m      <- lm(y_t[valid] ~ beta_vec[valid] + beta2_vec[valid] + sigma2_vec[valid])
  m_rob  <- coeftest(m, vcov = vcovHC(m, type = "HC3"))
  
  gammas_fm$gamma0[t] <- m_rob[1, 1]
  gammas_fm$gamma1[t] <- m_rob[2, 1]
  gammas_fm$gamma2[t] <- m_rob[3, 1]
  gammas_fm$gamma3[t] <- m_rob[4, 1]
}

gammas_fm <- gammas_fm %>% filter(!is.na(gamma0))
T_fm      <- nrow(gammas_fm)

# Fama-MacBeth average estimates and t-tests
fm_test <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  tv <- mu / se
  pv <- 2 * pt(-abs(tv), df = sum(!is.na(x)) - 1)
  c(mean = mu, se = se, t = tv, p = pv)
}

res <- sapply(list(
  gamma0 = gammas_fm$gamma0,
  gamma1 = gammas_fm$gamma1,
  gamma2 = gammas_fm$gamma2,
  gamma3 = gammas_fm$gamma3
), fm_test)

cat("\n=== Fama-MacBeth (1973) – Full Model Specification ===\n")
cat("Model: r_j = δ1 + δ2*β_j + δ3*β²_j + δ4*σ²(ε_j) + ε_j\n")
cat(sprintf("Number of monthly estimations: %d\n\n", T_fm))
cat(sprintf("δ1 (γ₀) = %8.5f  SE=%8.5f  t=%7.3f  p=%.4f\n",
            res["mean","gamma0"], res["se","gamma0"],
            res["t","gamma0"],    res["p","gamma0"]))
cat(sprintf("δ2 (γ₁) = %8.5f  SE=%8.5f  t=%7.3f  p=%.4f\n",
            res["mean","gamma1"], res["se","gamma1"],
            res["t","gamma1"],    res["p","gamma1"]))
cat(sprintf("δ3 (γ₂) = %8.5f  SE=%8.5f  t=%7.3f  p=%.4f\n",
            res["mean","gamma2"], res["se","gamma2"],
            res["t","gamma2"],    res["p","gamma2"]))
cat(sprintf("δ4 (γ₃) = %8.5f  SE=%8.5f  t=%7.3f  p=%.4f\n",
            res["mean","gamma3"], res["se","gamma3"],
            res["t","gamma3"],    res["p","gamma3"]))

cat("\nCAPM Validity Tests (δ3=0 and δ4=0):\n")
cat(sprintf("  δ3 (β²):    %s\n",
            ifelse(res["p","gamma2"] > 0.05,
                   "Fail to reject H0 (δ3=0) → Consistent with CAPM (Linearity)",
                   "Reject H0 (δ3=0) → Inconsistent with CAPM (Non-linearity)")))
cat(sprintf("  δ4 (σ²(ε)): %s\n",
            ifelse(res["p","gamma3"] > 0.05,
                   "Fail to reject H0 (δ4=0) → Consistent with CAPM (No ID risk pricing)",
                   "Reject H0 (δ4=0) → Inconsistent with CAPM (ID risk priced)")))

# Visualization of Gamma coefficients over time
df_gammas <- gammas_fm %>%
  mutate(date = dates_ret[t]) %>%
  pivot_longer(cols = c(gamma0, gamma1, gamma2, gamma3),
               names_to = "coefficient", values_to = "value") %>%
  mutate(coefficient = dplyr::recode(coefficient,
                                    gamma0 = "γ₀ (Intercept)",
                                    gamma1 = "γ₁ (Beta)",
                                    gamma2 = "γ₂ (Beta²) – Linearity Test",
                                    gamma3 = "γ₃ (σ²) – Non-systematic Risk Test"))

p_fm <- ggplot(df_gammas, aes(x = date, y = value)) +
  geom_line(color = "#003087", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#CC0000") +
  facet_wrap(~coefficient, scales = "free_y", ncol = 2) +
  labs(title = "Fama-MacBeth (1973) – Time Series of Gamma Coefficients",
       subtitle = "Red dashed line = zero (CAPM null hypothesis for γ₂ and γ₃)",
       x = NULL, y = "Coefficient Value") +
  theme_minimal(base_size = 9)

print(p_fm)
ggsave("fama_macbeth_coefficients.png", p_fm, width = 10, height = 6, dpi = 300)

# 5. DATA STORAGE -------------------------------------------------------------

save(capm_results, returns_list, mkt_prem, rf_ret,
     dates_ret, n100_vec, dates_index, rf_vec,
     portfolio_stats,
     file = "capm_data.RData")
