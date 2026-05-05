# Econometric Analysis of the CAPM on the Euronext 100 Index

This project provides a comprehensive econometric evaluation of the Capital Asset Pricing Model (CAPM) using a sample of 100 stocks from the Euronext 100 index between 2019 and 2025. The analysis covers beta estimation, cross-sectional testing of the Security Market Line (SML), Fama-MacBeth regressions, and structural stability, concluding with an asymmetric evaluation of market conditions.

## Project Structure

The project is executed via a series of R scripts, designed to be run sequentially using the provided master script.

*   `00_master_script.R`: The main script that runs the entire analysis sequentially.
*   `01_data_and_beta_estimation.R`: Handles data acquisition (Yahoo Finance, FRED), calculation of logarithmic returns, and Ordinary Least Squares (OLS) estimation of individual asset betas.
*   `02_capm_validity_testing.R`: Forms 10 beta-sorted portfolios, estimates the cross-sectional Security Market Line (SML), and performs the Fama-MacBeth (1973) two-pass regression to test for non-linearity and non-systematic risk pricing.
*   `03_structural_breaks_and_asymetry.R`: Investigates the temporal stability of beta using Chow tests and explores asymmetric beta behavior in up-markets versus down-markets.
*   `04_model_diagnostics.R`: Performs diagnostic checks on the underlying models (Breusch-Pagan, Shapiro-Wilk, Durbin-Watson, and Ramsey RESET).

## Data Sources

*   **Stock Prices and Market Index:** Monthly adjusted closing prices for 100 constituents of the Euronext 100 index and the index itself (`^N100`) sourced from Yahoo Finance.
*   **Risk-Free Rate:** The ECB Deposit Facility Rate (ECBDFR), sourced from the FRED database and converted to a monthly frequency.
*   **Time Horizon:** January 2019 to December 2025 (83 monthly observations).

## Key Findings

### 1. Beta Estimation and Jensen's Alpha

Individual regressions of the characteristic lines revealed that 94% of the assets had an alpha ( $\alpha$ ) that was not statistically significantly different from zero, which is largely consistent with the CAPM's assumption that market returns explain asset returns.

### 2. The Security Market Line (SML)

The cross-sectional regression of 10 portfolios formed on beta yielded the following SML:

$$\overline{r_p} - r_f = \gamma_0 + \gamma_1 \beta_p + \epsilon_p$$

The intercept ( $\gamma_0$ ) was found to be positive and statistically significant ( $p = 0.035$ ), contradicting the strict CAPM assumption that $\gamma_0 = 0$. The slope ( $\gamma_1$ ) was not significantly different from the realized market premium. This result — defensive stocks earning more than predicted and aggressive stocks earning less — is consistent with broader empirical literature.

### 3. Fama-MacBeth (1973) Regressions

The Fama-MacBeth methodology was applied to test for the linearity of the beta-return relationship and the pricing of non-systematic risk:

$$R_{pt} = \delta_{1t} + \delta_{2t} \beta_p + \delta_{3t} \beta_p^2 + \delta_{4t} \sigma^2(\epsilon_p) + \eta_{pt}$$

*   **Linearity:** The hypothesis $\delta_3 = 0$ could not be rejected, confirming a linear relationship between beta and returns.
*   **Non-Systematic Risk:** The hypothesis $\delta_4 = 0$ could not be rejected, confirming that the market does not reward idiosyncratic risk. Both findings support the CAPM framework.

### 4. Structural Breaks and Asymmetric Beta

*   **Stability:** Chow tests on a split sample (2019–2023 vs. 2024–2025) revealed that only 8.2% of the stocks exhibited a significant structural break in their parameters, suggesting beta is generally a stable parameter over time.
*   **Asymmetry:** An asymmetric CAPM specification differentiating between up-markets ( $D = 1$ ) and down-markets ( $D = 0$ ) was tested:

$$r_j - r_f = \alpha + \beta^{+} \cdot D \cdot (r_m - r_f) + \beta^{-} \cdot (1 - D) \cdot (r_m - r_f) + \epsilon$$

The asymmetric model showed substantially higher explanatory power ( $R^2 \approx 0.90$ ) compared to the symmetric model ( $R^2 = 0.43$ ). The slopes for both the up-market and down-market SMLs were highly significant, demonstrating that adjusting for market direction dramatically improves the model's fit.

## Diagnostics and Robustness

The analysis incorporated Newey-West and White HC3 robust standard errors to correct for heteroskedasticity and autocorrelation, ensuring reliable inference across the cross-sectional and time-series regressions.

## How to Run

1.  Clone this repository or download the `.R` scripts.
2.  Set your working directory in R/RStudio to the folder containing the scripts.
3.  Ensure the required packages are installed: `quantmod`, `tidyverse`, `ggplot2`, `gridExtra`, `strucchange`, `sandwich`, `lmtest`, `car`.
4.  Run the `00_master_script.R` file to execute the analysis pipeline sequentially.
