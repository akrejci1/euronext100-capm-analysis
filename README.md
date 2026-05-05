# CAPM Analysis — Euronext 100 (2019–2025)

This repository contains an empirical analysis of the Capital Asset Pricing Model (CAPM) on the Euronext 100 index over the period January 2019 – December 2025, carried out as part of an Econometrics course assignment.

---

## Overview

The project evaluates whether CAPM provides a satisfactory description of the cross-section of stock returns in the Euronext 100, with a focus on:

- Estimating beta coefficients and testing Jensen's alpha
- Validity of the Security Market Line (SML)
- Fama and MacBeth (1973) specification
- Stability of parameters over time (Chow tests)
- Beta asymmetry depending on market direction

---

## Data

| Variable | Source | Description |
|---|---|---|
| Stock prices | Yahoo Finance | 100 tickers, monthly closing prices |
| Market portfolio | Yahoo Finance (`^N100`) | Euronext 100 index |
| Risk-free rate | FRED (`ECBDFR`) | ECB deposit facility rate, converted to monthly frequency |
| Period | January 2019 – December 2025 | 83 monthly return observations |

Five tickers have incomplete data due to late listings or corporate events: `DSFIR.AS`, `URW.PA`, `EXO.AS`, `UMG.AS`, `JDEP.AS`.

---

## Repository Structure

```
.
├── 00_master_script.R                    # Runs all examples sequentially
├── 01_data_and_beta_estimation.R         # Data, beta estimation, alpha test
├── 02_capm_validity_testing.R            # SML, Fama-MacBeth test
├── 03_structural_breaks_and_asymmetry.R  # Chow tests, asymmetric beta
├── 04_model_diagnostics.R                # Diagnostics: BP, Shapiro-Wilk, DW, RESET, VIF
├── capm_beta_results.csv                 # Estimated beta coefficients for all 100 tickers
└── Results/                              # Output plots (PNG)
```

---

## Usage

### Requirements

```r
install.packages(c("quantmod", "tidyverse", "ggplot2",
                   "gridExtra", "strucchange", "lmtest",
                   "sandwich", "car"))
```

### Running the Analysis

Open `00_master_script.R` in RStudio and run the entire file. Scripts execute sequentially and pass data between each other via `capm_data.RData`.

Individual examples can also be run standalone — simply ensure `capm_data.RData` from the previous step is present in the working directory.

---

## Results

### 1 — Beta Estimation and Alpha Test

- Beta coefficients range from **0.15** (TEL.OL — Telenor) to **2.06** (A5G.IR — AIB Group)
- **41 defensive** stocks (β < 1), **59 aggressive** stocks (β ≥ 1)
- H₀: α = 0 rejected at 5% for **6 out of 100 tickers** (6.0%): KOG.OL, UNI.MI, PST.MI, BAMI.MI, SU.PA, ASM.AS — all with positive alpha

### 2 — CAPM Validity

**Security Market Line (OLS with White robust SE):**

| Parameter | Estimate | p-value | Interpretation |
|---|---|---|---|
| γ₀ (intercept) | 0.00517 | 0.035 | Rejected — inconsistent with CAPM |
| γ₁ (slope) | 0.00376 | 0.085 | Not rejected |
| Market premium | 0.00590 | — | Theoretical reference |

R² = 0.43 — beta explains 43% of cross-sectional variation in portfolio premia.

**Fama-MacBeth (1973) — monthly cross-sectional regressions (83 months):**

| Parameter | Estimate | p-value | CAPM hypothesis |
|---|---|---|---|
| δ₃ (β²) | −0.00146 | 0.703 | δ₃ = 0 ✓ — linearity confirmed |
| δ₄ (σ²) | 0.81073 | 0.287 | δ₄ = 0 ✓ — unsystematic risk not priced |

### 3 — Structural Breaks and Asymmetry

**Chow tests** (breakpoint: January 2024, 97 tickers analyzed):

- **8 tickers (8.2%)** show a significant structural break in both tests: UNI.MI, MC.PA, DSY.PA, STMPA.PA, SGO.PA, PST.MI, RI.PA, STLAM.MI
- **91.8%** of tickers show stable beta over time
- Both the structural break test and the predictive test yield identical results in all cases

**Asymmetric beta** (β⁺ for up-markets, β⁻ for down-markets, 50 up / 33 down months):

- **16 tickers (16.0%)** show statistically significant asymmetry at 5%
- Most asymmetric stocks have β⁻ > β⁺ (stronger reaction to market declines): AIR.PA, BIRG.IR, ENGI.PA
- Notable exceptions with β⁺ > β⁻: RACE.MI (Ferrari), MC.PA (LVMH)

**SML for up- and down-markets (Newey-West robust SE):**

| Parameter | Up-market | Down-market |
|---|---|---|
| γ₀ (p-value) | 0.00630 (p = 0.079) | −0.00694 (p = 0.041) |
| γ₁ (p-value) | **0.03091 (p < 0.001)** | **−0.02062 (p < 0.001)** |
| Market premium | 0.03277 | −0.03481 |
| R² | **0.930** | **0.876** |

The asymmetric CAPM dramatically outperforms the symmetric model (R² ≈ 0.90 vs. 0.43).

---

## Model Diagnostics

| Test | Statistic | p-value | Result |
|---|---|---|---|
| Breusch-Pagan (homoskedasticity) | 0.240 | 0.624 | ✓ |
| Shapiro-Wilk (normality) | 0.949 | 0.661 | ✓ |
| Durbin-Watson (no autocorrelation) | 1.869 | 0.263 | ✓ |
| Ramsey RESET (linearity) | 4.892 | 0.055 | ✓ |
| VIF — σ²(ε) in FM specification | 2.1 | — | ✓ |

High VIF for β and β² (≈ 22–26) is expected in the FM specification and does not affect the key coefficient on unsystematic risk.

---
