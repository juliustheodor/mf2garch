function [sum_alpha_beta] = alpha_beta_n_j(alpha, gamma, beta, j)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

sum_alpha_beta = 0;
		
for k = 1:j-2

sum_alpha_beta = sum_alpha_beta + ((alpha + gamma/2)+beta)^(k);

end
