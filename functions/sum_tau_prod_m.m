function [sum_tau_prod_mm] = sum_tau_prod_m(tau_forecast, t, s, alpha, gamma, beta, m)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

sum_tau_prod_mm = 0;

for j = 2:m 
	
    [sum_alpha_beta_j] = sum_alpha_beta(j, alpha, gamma, beta);
	sum_tau_prod_mm = sum_tau_prod_mm + tau_forecast(t+s-j) * sum_alpha_beta_j;
	
end

 

