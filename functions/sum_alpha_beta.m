function [sum_alpha_beta_j] = sum_alpha_beta(j, alpha, gamma, beta)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

sum_alpha_beta_j = 1; 

for k = 1:j-2
	
sum_alpha_beta_j = sum_alpha_beta_j + ((alpha + gamma/2)+beta)^k;
	
	
end