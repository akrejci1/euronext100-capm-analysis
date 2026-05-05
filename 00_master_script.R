# --- MASTER SCRIPT: CAPM ANALYSIS (EURONEXT 100) ---
# This script executes the full econometric pipeline sequentially.

# Clear environment, console, and graphics
rm(list = ls())
cat("\014")
graphics.off()

# Step 1: Data Acquisition and Beta Estimation
# Downloads data from Yahoo Finance/FRED and estimates individual asset betas.
source("01_data_and_beta_estimation.R", encoding = "UTF-8")

# Step 2: Testing CAPM Validity
# Constructs portfolios, estimates SML, and performs Fama-MacBeth regressions.
source("02_capm_validity_testing.R", encoding = "UTF-8")

# Step 3: Structural Breaks and Market Asymmetry
# Performs Chow tests for stability and analyzes up/down market beta asymmetry.
source("03_structural_breaks_and_asymmetry.R", encoding = "UTF-8")

# Step 4: Model Diagnostics
# Executes statistical tests for homoskedasticity, normality, and linearity.
source("04_model_diagnostics.R", encoding = "UTF-8")

cat("\n>>> ANALYSIS COMPLETE <<<\n")
cat("Final processed data saved in: capm_data.RData\n")
cat("All visualization outputs (PNG) are saved in the working directory.\n")
