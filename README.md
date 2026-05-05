# CAPM Analysis for Euronext 100 (2019–2025)

This repository contains an empirical analysis of the Capital Asset Pricing Model (CAPM) on the Euronext 100 index over the period 2019–2025.

## Overview

The project evaluates whether CAPM provides a satisfactory description of the cross‑section of stock returns in the Euronext 100, with a focus on beta estimates, the Security Market Line (SML), and potential asymmetries in risk–return relationships.  
The analysis is structured into three main parts: estimation of beta coefficients for individual stocks, testing the validity of CAPM using portfolio‑level regressions, and examining parameter stability and asymmetric betas over time.

## Data

- Monthly data from January 2019 to December 2025 are used, giving 83 monthly observations of logarithmic returns.  
- The investment universe consists of 100 stocks from the Euronext 100 index (Amsterdam, Brussels, Dublin, Lisbon, Milan, Oslo, Paris).  
- The market portfolio is proxied by the $$\(^\text{N100}\)$$ index.  
- The risk‑free rate is the ECB deposit facility rate (ECBDFR), converted to monthly frequency by dividing by 12.  
- Missing observations due to IPOs, mergers, and restructurings are handled by aligning time series with the market index; pre‑listing periods are excluded from regressions as NA.

## Methods

### 1. Individual stock CAPM regressions

For each stock $$\(j\)$$, the CAPM characteristic line is estimated using OLS on monthly excess returns:

$$\[r_j - r_f = \alpha_j + \beta_j (r_m - r_f) + \varepsilon_j\]$$

- At least 30 valid observations are required for each stock, which all 100 stocks satisfy.  
- Estimated betas range from approximately 0.15 (defensive, e.g. TEL.OL – Telenor) to 2.06 (aggressive, e.g. A5G.IR – AIB Group).  
- 41 stocks are classified as defensive $$(\(\beta < 1\))$$ and 59 as aggressive $$(\(\beta \ge 1\))$$.

Significance of abnormal returns is tested via:

- $$\(H_0: \alpha_j = 0\) vs. \(H_1: \alpha_j \neq 0\)$$ at the 5% level.  
- $$\(H_0\)$ is not rejected for 94 stocks (94%), while it is rejected for 6 stocks (all with positive alpha), a rate consistent with the expected number of false rejections at the 5% level.

### 2. Portfolio construction and SML

- The 100 stocks are sorted by estimated beta and grouped into 10 equally weighted portfolios (10 stocks each).  
- Portfolio beta is the arithmetic mean of constituent betas; portfolio excess return is the average monthly risk premium over the full sample.  

A cross‑sectional regression of portfolio average excess returns on portfolio betas is used to estimate the Security Market Line:

$$\[r_p - r_f = \gamma_0 + \gamma_1 \beta_p + \varepsilon_p\]$$

- Newey–West robust standard errors are applied due to the small sample of 10 portfolios.  
- The intercept $$\(\gamma_0\)$$ is positive and statistically significant (p ≈ 0.035), which is inconsistent with the CAPM restriction $$\(\gamma_0 = 0\)$$.  
- The slope $$\(\gamma_1\)$$ is not significant at the 5% level (p ≈ 0.085), but is close to the observed market risk premium, so $$\(\gamma_1 = r_m - r_f\)$$ cannot be rejected.  
- The regression explains about 43% of the cross‑sectional variation in portfolio risk premia $$(\(R^2 \approx 0.43\))$$.

These findings are in line with empirical literature: a positive, significant intercept and a slope lower than the market premium suggest that defensive stocks earn higher returns and aggressive stocks lower returns than predicted by CAPM.

### 3. Fama–MacBeth (1973) tests

The Fama–MacBeth procedure is implemented using cross‑sectional regressions across the 10 portfolios for each month $$\(t\)$$:

$$\[R_{pt} = \gamma_{0t} + \gamma_{1t} \beta_p + \gamma_{2t} \beta_p^2 + \gamma_{3t} \sigma^2(\varepsilon_p) + \eta_{pt}\]$$

- Time‑series averages of the monthly coefficients are computed, with standard errors robust to serial correlation.  
- Estimates are available for all 83 months.  
- The null of linearity in beta $$(\(\delta_3 = 0\)$$ for the $$\(\beta^2\)$$ term) is not rejected (p ≈ 0.703), supporting a linear beta–return relation.  
- The null that idiosyncratic risk is not priced $$(\(\delta_4 = 0\) for \(\sigma^2(\varepsilon)\))$$ is also not rejected (p ≈ 0.287).  

Multicollinearity between $$\(\beta\)$$ and $$\(\beta^2\)$$ is high (VIF ≈ 22–26), but this is expected and does not materially affect the key coefficient on idiosyncratic variance (VIF ≈ 2.1).

### 4. Parameter stability and structural breaks

To examine stability of CAPM parameters over time, the sample is split into two subperiods:

- Period 1: February 2019 – December 2023 (59 months).  
- Period 2: January 2024 – December 2025 (24 months).

A minimum of 20 valid observations is required in each period; three stocks (DSFIR.AS, URW.PA, EXO.AS) are excluded due to insufficient data, leaving 97 stocks.

Chow tests at the 5% level are applied for each stock:

- Structural‑break test: checks whether $$\(\alpha\)$$ and $$\(\beta\)$$ change between the two periods.  
- Forecast test: checks whether the model estimated in Period 1 predicts Period 2 returns well.

- 8 stocks (≈8.2%) show a significant structural break in both tests; 89 stocks (≈91.8%) do not.  
- Notable changes include LVMH (beta from ≈1.07 to ≈2.08) and Saint‑Gobain (beta from ≈1.63 to ≈0.36).  
- The high share of stocks without breaks indicates that beta is generally a stable parameter over this horizon.

SMLs are also estimated separately for each period:

- Period 1: intercept highly significant, slope not significant.  
- Period 2: intercept becomes insignificant (closer to CAPM prediction), while the slope is marginally significant, though the short sample (24 months, 10 portfolios) limits statistical power.

### 5. Asymmetric beta and market regimes

An asymmetric CAPM specification is used to allow different betas in up and down markets:

\[
r_j - r_f = \alpha + \beta^{+} D (r_m - r_f) + \beta^{-} (1 - D)(r_m - r_f) + \varepsilon
\]

where \(D = 1\) if the market rises and \(D = 0\) if it falls.

- The sample includes 50 months with rising markets (≈60.2%) and 33 months with falling markets (≈39.8%).  
- An F‑test at 5% significance is used to test \(\beta^{+} = \beta^{-}\).  
- 16 stocks (16%) exhibit statistically significant asymmetry; the remaining 84% have symmetric betas.  
- Most asymmetric stocks have \(\beta^{-} > \beta^{+}\), indicating higher sensitivity to downturns than to upturns (e.g. Airbus, Bank of Ireland, Engie), with Ferrari and LVMH as notable exceptions (\(\beta^{+} > \beta^{-}\)).

SMLs are then estimated separately for up and down markets using asymmetric betas and Newey–West standard errors:

- In rising markets: intercept ≈ 0.0063 (p ≈ 0.079), slope ≈ 0.0309 (p < 0.001), market premium ≈ 0.0328, \(R^2 \approx 0.93\).  
- In falling markets: intercept ≈ −0.0069 (p ≈ 0.041), slope ≈ −0.0206 (p < 0.001), market premium ≈ −0.0348, \(R^2 \approx 0.88\).  

The asymmetric CAPM dramatically improves fit relative to the symmetric model (\(R^2 \approx 0.90\) vs. 0.43), indicating that allowing for different betas in up and down markets captures much more of the cross‑sectional variation in returns.

## Main findings

- The majority of stocks have betas consistent with CAPM, and only a small fraction shows statistically significant abnormal returns; overall, results are broadly consistent with CAPM at the individual‑stock level.  
- The SML intercept is positive and significant, while the slope is close to but statistically indistinguishable from the market premium, suggesting deviations from the strict CAPM in the cross‑section of returns.  
- Fama–MacBeth tests support linearity in beta and confirm that idiosyncratic risk is not priced, in line with CAPM.  
- Betas are generally stable over time, with only about 8% of stocks exhibiting structural breaks between pre‑2024 and post‑2024 periods.  
- Roughly 16% of stocks display asymmetric betas, typically being more sensitive to market downturns than upturns.  
- The asymmetric CAPM specification provides a substantially better fit to the data than the symmetric version, suggesting that extensions of CAPM that account for asymmetry can better describe European equity returns in this sample.

## How to use this repository

- Clone the repository to reproduce the CAPM analysis for the Euronext 100 index.  
- Inspect the scripts and notebooks to see data preparation, model estimation, and test implementation (e.g. individual CAPM regressions, portfolio formation, SML, Fama–MacBeth, Chow tests, and asymmetric CAPM).  
- Adapt the code to different markets, time periods, or alternative risk factors if you want to extend the analysis beyond the basic and asymmetric CAPM.

## References

- Fama, E. F., & MacBeth, J. D. (1973). Risk, return, and equilibrium: Empirical tests.  
- Rachev, S. T. et al. (2007). Empirical evidence on the performance of CAPM and related models.
