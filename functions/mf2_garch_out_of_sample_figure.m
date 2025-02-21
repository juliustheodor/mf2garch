function  [] = mf2_garch_out_of_sample_figure(sigma_annual, an_vola_forecast, tau_forecast_annual, annual_unconditional_vola, foptions); 

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

% INPUTS: 
%   sigma_annual: (length(y) - 2*252 x 1 ) vector fitted values for annualized cond. volatility 
%   an_vola_forecast: (S x 1) vector of annualized cond. volatility forecasts 
%   tau_forecast_annual: (length(y) - 2*252 + S x 1 ) vector fitted values for
%       annualized long-term component extend using forecasts 
%   annual_unconditional_vola: scalar of annualized unconditional volatility 
%   foptions

% OUTPUT: 
% Figure saved as 'ForecastEndofSample.png' in the figures folder.

% Construct variable used for x-axis of figure: 
Timing = (-50:foptions.S)';

% Figure for out-of-sample forecasting: 
figure
plot(Timing(1:50),sigma_annual(end-49:end),'k-','LineWidth',3); hold on;
plot(Timing(51:50+foptions.S), an_vola_forecast(1:foptions.S), 'k--','LineWidth',2);
plot(Timing(1:50), tau_forecast_annual(end-49-foptions.S:end-foptions.S),'r-', 'LineWidth',2);
plot(Timing(51:51+foptions.S), tau_forecast_annual(end- foptions.S :end),'r--','LineWidth',2);
yline(annual_unconditional_vola, 'b--'); 
xlabel('Date', 'Interpreter', 'latex'); ylabel('conditional/long-term volatility', 'Interpreter', 'latex'); hold off 
xline(0, '--', 'LineWidth', 1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
xlim([-50 foptions.S]);
legend({'$\sigma_t$ conditional volatility (in-sample)', ...
        '$\sigma_t$ conditional volatility (out-of-sample forecast)', ... 
        '$\tau_t$ long-term volatility (in-sample)', ...
        '$\tau_t$ long-term volatility (out-of-sample forecast)', ...
        'unconditional volatility'
        }, ...
        'Interpreter', 'latex', 'Location', 'best');
set(gcf, 'Color', 'w');
set(gca,'TickLabelInterpreter','latex', 'FontSize', 14); 
saveas(gcf,'figures/ForecastEndofSample.png'); 

end 