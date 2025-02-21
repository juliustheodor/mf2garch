function [delta_m] = Delta_m(alpha, gamma, beta, lambda_0, lambda_1, lambda_2, kappa, m)
% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

sum_1 = 0;
	
for j = 2:m
	
[sum_alpha_beta] = alpha_beta_n_j(alpha, gamma, beta, j);    
    
sum_1 = sum_1 + sum_alpha_beta;
	
end
    
delta_m = (1-(alpha + gamma/2)-beta)* lambda_1 * ((alpha + gamma/2)+beta) * lambda_0/(1-lambda_1-lambda_2) * ( (m-1)/m + 1/m * sum_1);

