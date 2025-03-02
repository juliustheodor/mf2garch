%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          MF2 GARCH-rw-m Model
%                   Christian Conrad & Julius Schoelkopf 
%                   Heidelberg University, February 2025 
%                               Version 0.1.0
%
% This file provides code for four applications for estimating the 
% MF2-GARCH-rw-m model and computing volatility forecasts.
% 
% Example A: Estimation of the MF2-GARCH-rw-m using S&P500 return data
%            including a plot of the time series of fitted conditional
%            variances 
% Example B: News impact curve (NIC) 
% Example C: Out-of-sample forecasting
% Example D: Illustration of forecast behavior 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

addpath 'data/'
addpath 'functions/'
addpath 'figures/'

%% Import the return data to Matlab (S&P500 returns from 1971-2023) 

% Read the data into a table

Returns = readtable('data/SP500_1971_2023_06_30_ret.xlsx');

% Extract the column 'RET_SPX' from the table and store it 

y = Returns.RET_SPX;

%% Select lag length m 
% For the long-term component, specify m, i.e. the number days over which
% V_{t-1}^m is computed. Choose whether you want to use a fixed value of m or 
% let the optimal m be selected as the one that minimizes the BIC:         

foptions.choice = 'BIC'; % choices: 'BIC' or 'fix' (specify m) 

% If you select `BIC', the code will save the optimal m after running 
% mf2_garch_estimation for forecasting or the NIC. Moreover, the code will
% generate a figure that plots the BIC as a function of m. 

% If f.options.choice = 'fix', please specify the m you choose here: 

foptions.m=63; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Example A (Estimation and Time Series) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Estimation of the MF2-GARCH-rw-m 
% The following line of code replicates the second panel of Table 2 in 
% Conrad & Engle (2025, doi.org/10.1002/jae.3118) for the MF2-GARCH-rw-m. 
% In Conrad & Engle (2025), all models were estimated using OxMetrics. 

% The function gives you an estimation output in the command window as 
% well as vectors for the coefficients, the standard errors, the p-values, 
% the standardized residuals, the fitted values for the short and long-term component 
% as well as the time series for the (annualized) conditional variance and 
% the estimate for the annualized unconditional volatility. 

% Bollerslev-Wooldridge robust standard errors are reported.
% The function uses constraints on the parameters following Assumption 2 (for the 
% short-term component) and Assumption 3 (for the long-term component) of 
% Conrad & Engle (2025). When computing the likelihood, we discard the first two 
% years of the return series (i.e., 2*252 trading days) to account for lags 
% of the squared deGARCHed returns in the long-term component. This allows 
% comparing the BIC of models with different values of m. 
% For details on the estimation, see section A.1.1 in Conrad & Engle (2025). 

[coeff, qmle_se, p_value_qmle,  Z, h, tau, sigma_annual, tau_annual, annual_unconditional_vola, foptions]  = mf2_garch_estimation(y,foptions); 

% The default here is the full return sample. You can also use a subsample of y. 

%% Plot of time series 
% Extract the date column (not required for estimation, only for figure) 

dates = datetime(Returns.OBS, 'InputFormat', 'MM/dd/yyyy'); 

% The figure shows the estimated conditional volatility and long-term volatility
% over the full-sample. All quantities are annualized.
% Grey shaded areas represent NBER recession periods in the US.

mf2_garch_time_series(dates, sigma_annual, tau_annual); 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Example B (News impact curve) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the news impact curve for the estimated model 
% The NIC is presented in terms of annualized volatilities (see equation (10) 
% in Conrad & Engle (2025)). 

[r, NIC] = mf2_garch_nic(Z, h, tau, foptions, coeff);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Example C (Forecasting out of sample) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specifiy the maximum forecasting horizon: 

foptions.S = 252;

% Forecasting at the end of the sample S (trading) days into the future: 
% You must use the same sample as in estimation function!

[horizon, forecast, an_vola_forecast, h_forecast, tau_forecast, tau_forecast_annual]  = mf2_garch_forecasting(y, Z, h, tau, coeff, foptions);

% The function displays the forecasts (from the end of the sample) for the 
% annualized volatility on the next day, next week (5 days), next month 
% (21 days), next 6 months (126 days), and next 12 months (252 days) based 
% on the estimated parameters.

% The following code yields a figure of the forecasts of long-term 
% volatility and conditional volatility. 
% You must use the same sample as in estimation function!

mf2_garch_out_of_sample_figure(sigma_annual, an_vola_forecast, tau_forecast_annual, annual_unconditional_vola, foptions)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Example D (Illustration of forecast behaviour) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete previous fitted values and coefficients because we use a new
% estimation sample: 

clear tau h coeff qmlse_se  Z  tau_annual annual_unconditional_vola

% Specifiy the maximum forecasting horizon: 

foptions.S = 120; 

% Estimation of the MF2-GARCH. We want to forecast volatility from August 10, 
% 2011 (10249 in dates vector) S = 120 days into the future. 
% Specify the cutoff from where you want to forecast: 

foptions.cutoff_date = datetime(2011,8,10);  
foptions.cutoff = 10249; 

% Therefore, we need to reestimate the model using data until August 10, 2011. 

[coeff, ~, ~, Z, h, tau, ~, tau_annual, annual_unconditional_vola, foptions]  = mf2_garch_estimation(y(1:foptions.cutoff),foptions); 

% Forecasting exercise 
% This function provides forecasts for the annualized volatility, h and tau 
% for the next S days from the end of the specified sample. 

[horizon, forecast, an_vola_forecast, h_forecast, tau_forecast, tau_forecast_annual]  = mf2_garch_forecasting(y(1:foptions.cutoff), Z, h, tau, coeff, foptions);

% Illustration of forecasting behaviour as in Figure 5 in Conrad & Engle (2025): 

mf2_garch_illustration_forecasting_figure(sigma_annual, an_vola_forecast, tau_forecast_annual, annual_unconditional_vola, foptions, dates)