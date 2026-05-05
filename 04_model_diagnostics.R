# --- STEP 4: MODEL DIAGNOSTICS ---
# This script performs statistical tests to verify CAPM model assumptions.
# Author: Adam Krejci

library(lmtest)
library(car)

# Ensure data is loaded
if (!exists("sml_model")) load("capm_data.RData")

cat("\n=== 1. MAIN SML MODEL DIAGNOSTICS (FULL PERIOD) ===\n")

# 1. HETEROSKEDASTICITY TEST (Breusch-Pagan)
# H0: Constant variance of errors (homoskedasticity)
bp_test <- bptest(sml_model)
cat(sprintf("Test 1: Breusch-Pagan test: BP = %.3f, p-value = %.4f\n", 
            bp_test$statistic, bp_test$p.value))
cat("Interpretation:", ifelse(bp_test$p.value > 0.05, 
                            "Status: OK (Homoskedasticity)", 
                            "Warning: Heteroskedasticity detected (non-constant error variance)"), "\n")

# 2. NORMALITY OF RESIDUALS (Shapiro-Wilk)
# H0: Residuals are normally distributed
sw_test <- shapiro.test(residuals(sml_model))
cat(sprintf("Test 2: Shapiro-Wilk test: W = %.3f, p-value = %.4f\n", 
            sw_test$statistic, sw_test$p.value))
cat("Interpretation:", ifelse(sw_test$p.value > 0.05, 
                            "Status: OK (Normally distributed residuals)", 
                            "Warning: Non-normal residuals (common in financial time series)"), "\n")

# 3. AUTOCORRELATION TEST (Durbin-Watson)
# H0: Residuals are not autocorrelated (DW target ≈ 2)
dw_test <- dwtest(sml_model)
cat(sprintf("Test 3: Durbin-Watson test: DW = %.3f, p-value = %.4f\n", 
            dw_test$statistic, dw_test$p.value))
cat("Interpretation:", ifelse(dw_test$p.value > 0.05, 
                            "Status: OK (No autocorrelation)", 
                            "Warning: Autocorrelation detected in residuals"), "\n")

# 4. RAMSEY RESET TEST (Linearity)
# H0: The linear functional form is correct
reset_main <- resettest(sml_model, power = 2:3, type = "regressor")
cat(sprintf("Test 4: Ramsey RESET test: F = %.3f, p-value = %.4f\n", 
            reset_main$statistic, reset_main$p.value))
cat("Interpretation:", ifelse(reset_main$p.value > 0.05, 
                            "Status: OK (Relationship is linear)", 
                            "Warning: Possible non-linearity detected"), "\n")


cat("\n=== 2. SUB-MODEL DIAGNOSTICS (BREUSCH-PAGAN TESTS) ===\n")

# 5. Period 1 (Stable)
bp_sml1 <- bptest(sml_model_1)
cat(sprintf("Test 5: Period 1: BP = %.3f, p-value = %.4f\n", bp_sml1$statistic, bp_sml1$p.value))
cat("Interpretation:", ifelse(bp_sml1$p.value > 0.05, "Status: OK", "Warning: Heteroskedasticity detected"), "\n")

# 6. Period 2 (Last 24m)
bp_sml2 <- bptest(sml_model_2)
cat(sprintf("Test 6: Period 2: BP = %.3f, p-value = %.4f\n", bp_sml2$statistic, bp_sml2$p.value))
cat("Interpretation:", ifelse(bp_sml2$p.value > 0.05, "Status: OK", "Warning: Heteroskedasticity detected"), "\n")

# 7. Up-Market Condition
bp_up <- bptest(sml_m_up)
cat(sprintf("Test 7: Up-Market: BP = %.3f, p-value = %.4f\n", bp_up$statistic, bp_up$p.value))
cat("Interpretation:", ifelse(bp_up$p.value > 0.05, "Status: OK", "Warning: Heteroskedasticity detected"), "\n")

# 8. Down-Market Condition
bp_dn <- bptest(sml_m_dn)
cat(sprintf("Test 8: Down-Market: BP = %.3f, p-value = %.4f\n", bp_dn$statistic, bp_dn$p.value))
cat("Interpretation:", ifelse(bp_dn$p.value > 0.05, "Status: OK", "Warning: Heteroskedasticity detected"), "\n")


cat("\n=== 3. MULTICOLLINEARITY (VIF) FOR FAMA-MACBETH SPECIFICATION ===\n")

# 9. VIF Calculation
vif_model <- lm(avg_prem ~ beta_p + I(beta_p^2) + sigma2_eps, data = portfolio_stats)
vif_values <- vif(vif_model)
cat("Test 9: VIF Values:\n")
print(vif_values)
cat("Note: VIF > 10 for beta and beta^2 is common and expected in Fama-MacBeth models.\n")
cat("Focus: Low VIF for sigma2_eps (idiosyncratic risk) is the primary requirement.\n")

cat("\n--- DIAGNOSTICS COMPLETE ---\n")
