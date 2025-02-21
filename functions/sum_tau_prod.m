function [sum_tau_prod_s] = sum_tau_prod(tau_forecast, t, s, alpha, gamma, beta)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

sum_tau_prod_s = 0;

for j = 2:s-1
	
    [sum_alpha_beta_j] = sum_alpha_beta(j, alpha, gamma, beta);
	sum_tau_prod_s = sum_tau_prod_s + tau_forecast(t+s-j) * sum_alpha_beta_j;
	
end