[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# MF2-GARCH Toolbox for Matlab (developed by Christian Conrad and Julius Schoelkopf, 2025)

A Matlab package for estimating and forecasting using the multiplicative factor multi-frequency GARCH (MF2-GARCH)  as proposed in Conrad & Engle (2025) accompanying the paper „ Long-term volatility shapes the stock markets sensitivity to news“ by Conrad, Schoelkopf, and Tushteva (2024): 

* A comprehensive toolbox for estimating and forecasting using the MF2-GARCH-rw-m.
* Five applications: estimation, news-impact-curve, illustration of long-term component, out-of-sample forecasting, illustration of forecasting behavior 


## Installation



## Suggested Citation
Please cite as: 
> Conrad, Christian and Schoelkopf, Julius Theodor and Tushteva, Nikoleta, Long-Term Volatility Shapes the Stock Market's Sensitivity to News (2024). Available at SSRN:  http://dx.doi.org/10.2139/ssrn.4632733

and 

> Conrad, Christian and Julius Schoelkopf. 2025. MF2-GARCH Toolbox for Matlab. Matlab package version 0.1.0. 

## Contact 
Please address any questions about the Matlab code to:
Julius Schoelkopf, Heidelberg University, Department of Economics. Email: julius.schoelkopf [at] awi.uni-heidelberg.de 

We do not assume any responsibilities for results produced with the available code. Please let me know, if you have suggestions for further versions or find any bugs. 

# Example 

The following line of code replicates the second panel in Table 2 in Conrad & Engle (2025) for the MF2-GARCH-rw-m. In Conrad & Engle (2025) all models were estimates using OxMetrics. Bollerslev-Wooldridge robust standard errors are reported. The Matlab function uses constraints on the parameters following assumption 2 (for the short-term component) and assumption 3 (for the long-term component) of Conrad & Engle (2025). For the fitted values, we discard the first two years of y (i.e., 2 times 252 trading days) to account for lags of the squared deGARCHed returns when comparing models using the BIC. You could decrease this, but you need to discard at least 2m values. For details on the estimation, see section A.1.1 in Conrad & Engle (2025). 

```matlab
[coeff, qmle_se, p_value_qmle,  Z, h, tau, sigma_annual, tau_annual, annual_unconditional_vola, foptions]  = mf2_garch_estimation(y,foptions); 
```
