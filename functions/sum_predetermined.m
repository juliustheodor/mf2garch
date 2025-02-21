function [sum_predetermined_m] = sum_predetermined(r2, h_forecast, s, t, m)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

sum_predetermined_m = 0;

for j = s:m

sum_predetermined_m = sum_predetermined_m + r2(t+s-j)./h_forecast(t+s-j);

end