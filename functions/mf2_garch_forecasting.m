function [horizon, forecast, an_vola_forecast, h_forecast, tau_forecast, tau_forecast_annual] = mf2_garch_forecasting(y, Z, h, tau, coeff, foptions)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

% INPUTS: 
%   y: a (Tx1) vector with the return time series 
%   Z: (length(y) - 2*252 x 1 ) vector with (Z= (y-mu) ./ sqrt(h.*tau)) 
%   h: (length(y) - 2*252 x 1 ) vector of fitted values for short-term component 
%   tau: (length(y) - 2*252 x 1 ) vector fitted values for long-term component 
%   coeff: (7x1) coefficient vector for the MF2-GARCH, estimation results for
%       (mu, alpha, gamma, beta, lambda_0, lambda_1, lambda_2)' 
%   foptions: 
%       f.options.m: a scalar (used for the rolling window measure of the local bias of 
%       the short-term component’s conditional variance over the previous m days) 
%       For forecasting, you must specify the m you want to use using f.options.m. 
%       f.options.S: a scalar that determines the maximum forecasting horizon.  

% OUTPUTS: 
%   horizon: (length(y) + Sx1) vector with horizon counting 
%   forecast: length(y) - 2*252 + S x 1 ) vector of fitted values for cond. volatility
%   an_vola_forecast: (S x 1) vector of annualized cond. volatility forecasts 
%   h_forecast: (length(y) - 2*252 + S x 1 ) vector of fitted values for short-term component 
%   tau_forecast: (length(y) - 2*252 + S x 1 ) vector fitted values for long-term component 
%   tau_forecast_annual (length(y) - 2*252 + S x 1 ) vector fitted values
%       for long-term component (annualized) 


    % parameters: 
    
    mu = coeff(1);
    alpha = coeff(2);
    gamma = coeff(3);
    beta = coeff(4);
    
    lambda_0 = coeff(5);
    lambda_1 = coeff(6);
    lambda_2 = coeff(7);
    
    kappa = mean((Z.^2 - 1).^2) + 1; 

    % number of observations 
	t = length(Z); 
	
    m = foptions.m; 

    % maximum forecasting horizon 
    S = foptions.S; 

	% Compute unconditional variance
    
    [delta_m] = Delta_m(alpha, gamma, beta, lambda_0, lambda_1, lambda_2, kappa, m);

    % starting values for r2, h and tau: observations 1 to t
	
    leng = length(Z)+S;
    horizon  = zeros(leng,1);
	r2       = zeros(leng,1); 
	tau_forecast      = zeros(leng,1);
	h_forecast        = zeros(leng,1);
	forecast = ones(leng,1);
    
    r2(1:t)           = Z.^2 .* h .* tau;
	h_forecast(1:t)   = h;
	tau_forecast(1:t) = tau;
	forecast(1:t)  = h.*tau;
    
    % h and tau in t+1
	if (y(length(y)) >= 0 )
	
	h_forecast(t+1) = (1-(alpha + gamma/2)-beta) + alpha * r2(t)./tau_forecast(t) + beta * h_forecast(t) ;
    
    end
    
	if (y(length(y)) < 0 )
	
	h_forecast(t+1) = (1-(alpha + gamma/2)-beta) + (alpha + gamma )* r2(t)./tau_forecast(t) + beta * h_forecast(t) ;
    
    end
    
    [sum_predetermined_m] = sum_predetermined(r2, h_forecast, 1, t, m);
	tau_forecast(t+1) = lambda_0 + lambda_1 * 1/m * sum_predetermined_m + lambda_2 * tau_forecast(t);
    
    % Forecasting the GARCH component
	for k = t+1:t+S
	
	s=k-t;
	horizon(t+s) = s;
	h_forecast(t+s) = 1+((alpha + gamma/2)+beta)^(s-1)*(h_forecast(t+1)-1);
	
    end
    
    % Forecasting the long-term component	
	if (S<=m)
	
	for k = t+2:t+S
	
	s=k-t;
    [sum_tau_s_1] = sum_tau(tau_forecast, s, t);
    [sum_predetermined_m] = sum_predetermined(r2, h_forecast, s, t, m);
	tau_forecast(t+s) = lambda_0 + (lambda_1 * 1/m + lambda_2)*tau_forecast(t+s-1) + lambda_1 * 1/m * sum_tau_s_1 + lambda_1 * 1/m * sum_predetermined_m;
    
    end

    end

	if (S>m)
	
	for k = t+2:t+m
	
	s=k-t;
    [sum_tau_s_1] = sum_tau(tau_forecast, s, t);
    [sum_predetermined_m] = sum_predetermined(r2, h_forecast, s, t, m);
	tau_forecast(t+s) = lambda_0 + (lambda_1 * 1/m + lambda_2)*tau_forecast(t+s-1) + lambda_1 * 1/m * sum_tau_s_1 + lambda_1 * 1/m * sum_predetermined_m;
    
    end

	for k = t+m+1:t+S
	
	s=k-t;
    [sum_tau_all_m] = sum_tau_all(tau_forecast, s, t, m);
	tau_forecast(t+s) = lambda_0 + (lambda_1 * 1/m + lambda_2)*tau_forecast(t+s-1) + lambda_1 * 1/m * sum_tau_all_m;
    
    end
	
    end
    
    % Forecasting the conditional variance

	forecast(t+1) = h_forecast(t+1) .* tau_forecast(t+1);

	if (S<=m)
	
	for k = t+2:t+S
	
	s=k-t;
    [sum_predetermined_m] = sum_predetermined(r2, h_forecast, s, t, m);
    [sum_tau_prod_s] = sum_tau_prod(tau_forecast, t, s, alpha, gamma, beta);
    [sum_forecast_j] = sum_forecast(forecast, alpha, gamma, beta, t, s);
	forecast(t+s) = (1-(alpha + gamma/2)-beta) * tau_forecast(t+s) + lambda_0 * ((alpha + gamma/2)+beta) * h_forecast(t+s-1) +  (lambda_1 * 1/m * ((alpha + gamma/2) * kappa +beta) + lambda_2 * ((alpha + gamma/2)+beta)) * forecast(t+s-1) +  lambda_1 * ((alpha + gamma/2)+beta) * h_forecast(t+s-1) * 1/m * sum_predetermined_m + (1-(alpha + gamma/2)-beta) * lambda_1 * ((alpha + gamma/2)+beta) * 1/m * sum_tau_prod_s +  lambda_1 * ((alpha + gamma/2)*kappa+beta) * ((alpha + gamma/2)+beta) * 1/m * sum_forecast_j ;
    end

    end

 	if (S>m)
	
	for k = t+2:t+m 
	
	s=k-t;
    [sum_predetermined_m] = sum_predetermined(r2, h_forecast, s, t, m);
    [sum_tau_prod_s] = sum_tau_prod(tau_forecast, t, s, alpha, gamma, beta);
    [sum_forecast_j] = sum_forecast(forecast, alpha, gamma, beta, t, s);
	forecast(t+s) = (1-(alpha + gamma/2)-beta) * tau_forecast(t+s) + lambda_0 * ((alpha + gamma/2)+beta) * h_forecast(t+s-1) +  (lambda_1 * 1/m * ((alpha + gamma/2) * kappa +beta) + lambda_2 * ((alpha + gamma/2)+beta)) * forecast(t+s-1) +  lambda_1 * ((alpha + gamma/2)+beta) * h_forecast(t+s-1) * 1/m * sum_predetermined_m + (1-(alpha + gamma/2)-beta) * lambda_1 * ((alpha + gamma/2)+beta) * 1/m * sum_tau_prod_s +  lambda_1 * ((alpha + gamma/2)*kappa+beta) * ((alpha + gamma/2)+beta) * 1/m * sum_forecast_j ;
    end


	for k = t+m+1:t+S
	
	s=k-t;
    [sum_tau_prod_mm] = sum_tau_prod_m(tau_forecast, t, s, alpha, gamma, beta, m);
    [sum_forecast_mm] = sum_forecast_m(forecast, alpha, gamma, beta, t, s, m);
	forecast(t+s) = (1-(alpha + gamma/2)-beta) * tau_forecast(t+s) + lambda_0 * ((alpha + gamma/2)+beta) * h_forecast(t+s-1) +  (lambda_1 * 1/m * ((alpha + gamma/2) * kappa +beta) + lambda_2 * ((alpha + gamma/2)+beta)) * forecast(t+s-1) + (1-(alpha + gamma/2)-beta) * lambda_1 * ((alpha + gamma/2)+beta) * 1/m * sum_tau_prod_mm +  lambda_1 * ((alpha + gamma/2)*kappa+beta) * ((alpha + gamma/2)+beta) * 1/m * sum_forecast_mm ;
    end

    end

    % Annualize forecasts for tau 
    tau_forecast_annual = sqrt(252*tau_forecast);
    
    an_vola_forecast = (252 .* forecast(t+1:t+S)).^0.5;
    disp('annualized volatility forecast 1 day: ');
    disp(an_vola_forecast(1));
    if S > 4 
    disp('annualized volatility forecast 1 week (5 days): ');
    disp(an_vola_forecast(5));
    end 

    if S > 20 
    disp('annualized volatility forecast 1 moth (21 days): ');
    disp(an_vola_forecast(21));
    end 

    if S > 125 
    disp('annualized volatility forecast 6 moths (126 days): ');
    disp(an_vola_forecast(126));
    end 

    if S > 251 
    disp('annualized volatility forecast 1 year (252 days): ');
    disp(an_vola_forecast(252));
    end 
end 