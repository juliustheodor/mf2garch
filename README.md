# MF2-GARCH Toolbox for Matlab (developed by Christian Conrad and Julius Schoelkopf, 2025)

A Matlab package for estimating and forecasting using the multiplicative factor multi-frequency GARCH (MF2-GARCH)  as proposed in Conrad & Engle (2025) accompanying the paper „Long-term volatility shapes the stock markets sensitivity to news“ by Conrad, Schoelkopf, and Tushteva (2024): 

* A comprehensive toolbox for estimating and forecasting using the MF2-GARCH-rw-m.
* Five applications: estimation, news-impact-curve, illustration of long-term component, out-of-sample forecasting, illustration of forecasting behavior 

## Suggested Citation
Please cite as: 
> Conrad, Christian and Schoelkopf, Julius Theodor and Tushteva, Nikoleta, Long-Term Volatility Shapes the Stock Market's Sensitivity to News (2024). Available at SSRN:  http://dx.doi.org/10.2139/ssrn.4632733

and 

> Conrad, Christian and Julius Schoelkopf. 2025. MF2-GARCH Toolbox for Matlab. Matlab package version 0.1.0. 

## Contact 
Please address any questions about the Matlab code to:
* Christian Conrad,  Heidelberg University, Department of Economics. Email: christian.conrad [at] awi.uni-heidelberg.de 
* Julius Schoelkopf, Heidelberg University, Department of Economics. Email: julius.schoelkopf [at] awi.uni-heidelberg.de 

We do not assume any responsibilities for results produced with the available code. Please let us know, if you have suggestions for further versions or find any bugs. 

# Applications 

## Estimation of the MF2-GARCH-rw-m model in Matlab for S&P 500 stock returns. 
The following application of the MF-2GARCH replicates the second panel in Table 2 in Conrad & Engle (2025) for the MF2-GARCH-rw-m. In Conrad & Engle (2025), all models were estimates using OxMetrics. The Matlab function uses constraints on the parameters following assumption 2 (for the short-term component) and assumption 3 (for the long-term component) of Conrad & Engle (2025). 

Define daily log-returns as $y_t=\sigma_t Z_t=$ $\sqrt{h_t \tau_t} Z_t$ where $Z_t$ is i.i.d. and has a symmetric density with mean zero and variance one. $\sigma_t^2$ denotes the conditional variance and the short- and long-term volatility components are given by $h_t$ and $\tau_t$. Let `y` be this (Tx1) vector of daily log-returns. 

The short-term volatility component is defined as a unit variance GJR-GARCH(1,1)
$$
h_t=(1-\phi)+\left(\alpha+\gamma \mathbf{1}_{\left\{r_{t-1}<0\right\}}\right) \frac{r_{t-1}^2}{\tau_{t-1}}+\beta h_{t-1}
$$
and the long-term component is specified as a MEM equation for the conditional expectation of $V_t$ (squared deGARCHed returns):
$$
\tau_t=\lambda_0+\lambda_1 V_{t-1}^{(m)}+\lambda_2 \tau_{t-1}
$$
where
$$
V_{t-1}^{(m)}=\frac{1}{m} \sum_{j=1}^m V_{t-j}=\frac{1}{m} \sum_{j=1}^m \frac{r_{t-j}^2}{h_{t-j}}.
$$

The function `mf2_garch_estimation(y,foptions)` provides you with an estimation output for the seven parameters of the short- and long-term component in the command window as well as fitted values for the conditional volatility as well as the long- and short-term component. The output of the function are the vectors for the coefficient estimates (`coeff')`, the Bollerslev-Wooldridge robust standard errors  (`qmle_se`), the p-values  (`p_value_qmle`), the shocks (`Z`), the fitted values for the short (`h`) and long-term component (`tau` or annualized `tau_annual`) as well as the time series for the annualized conditional volatility (`sigma_annual`) and the estimate for the annualized unconditional volatility (`annual_unconditional_vola`). For the long-term component, you need to specify m, i.e. the number days over which $V_t^m$ is computed. Choose whether you want to use a fixed value of m or let the optimal m be selected as the one that minimizes the BIC. The `foptions` structure contains the researcher's choice for m. You either specifiy `foptions.choice = 'BIC'` if the optimal $m$ needs to be selected or `foptions.choice = 'fix'` together with the choice of your m as `foptions.m=63`. For the fitted values, the code discards the first two years of y (i.e., 2 times 252 trading days) to account for lags of the squared deGARCHed returns when comparing models using the BIC. You could decrease this, but you need to discard at least 2m values. For details on the estimation, see section A.1.1 in Conrad & Engle (2025). 

```matlab
[coeff, qmle_se, p_value_qmle,  Z, h, tau, sigma_annual, tau_annual, annual_unconditional_vola, foptions]  = mf2_garch_estimation(y,foptions); 
```

<img src="figures/ForecastEndofSample.png" width="80%" />
