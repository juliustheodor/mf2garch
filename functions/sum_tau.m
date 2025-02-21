function [sum_tau_s_1] = sum_tau(tau_forecast, s, t)

sum_tau_s_1 = 0;
	
for j = 2:s-1
	
sum_tau_s_1 = sum_tau_s_1 + tau_forecast(t+s-j);
	
end

