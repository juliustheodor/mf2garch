function [sum_forecast_mm] = sum_forecast_m(forecast, alpha, gamma, beta, t, s, m)

sum_forecast_mm = 0;

for j = 2:m
	
	sum_forecast_mm = sum_forecast_mm + ((alpha + gamma/2)+beta)^(j-2) * forecast(t+s-j);
	
end
