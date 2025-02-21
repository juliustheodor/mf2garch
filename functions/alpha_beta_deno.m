function [sum_alpha_beta_deno] = alpha_beta_deno(alpha, gamma, beta, m)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

sum_alpha_beta_deno = 0;
		
for j = 2:m

sum_alpha_beta_deno = sum_alpha_beta_deno + ((alpha + gamma/2)+beta)^(j-1);

end