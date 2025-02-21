function [sum_tau_all_m] = sum_tau_all(tau_forecast, s, t, m)

 sum_tau_all_m = 0;
	
for j = 2:m
	
	sum_tau_all_m = sum_tau_all_m + tau_forecast(t+s-j);

end