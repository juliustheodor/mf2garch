function [sum_forecast_j] = sum_forecast(forecast, alpha, gamma, beta, t, s)

sum_forecast_j = 0;

for j = 2:s-1

sum_forecast_j = sum_forecast_j + ((alpha + gamma/2)+beta)^(j-2) * forecast(t+s-j);

end
	
