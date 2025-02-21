function  [] = mf2_garch_illustration_forecasting_figure(sigma_annual, an_vola_forecast, tau_forecast_annual, annual_unconditional_vola, foptions, dates) 

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

% INPUTS: 
%   an_vola_forecast: (S x 1) vector of annualized cond. volatility forecasts 
%   sigma_annual: (length(y) - 2*252 x 1 ) vector fitted values for annualized cond. volatility 
%   tau_annual: (length(y) - 2*252 x 1 ) vector fitted values for annualized long-term component 
%   annual_unconditional_vola: scalar of annualized unconditional volatility 
%   foptions:
%       foptions.cutoff_date: date from which you want to forecast 
%       foptions.cutoff: number of row in date vector from which you want to forecast 
%   dates: vector of dates (not adjusted as tau and sigma, thus we use 2*252 in code below) 


% OUTPUT: 
% Figure saved as 'ForecastIllustration.png' in the figures folder.

cutoff_date = foptions.cutoff_date; 
cutoff = foptions.cutoff; 

figure
plot(dates(cutoff-49:cutoff+foptions.S),sigma_annual(cutoff-2*252-49:cutoff-2*252+foptions.S),'k-','LineWidth',3); hold on;
plot(dates(cutoff+1:cutoff+foptions.S), an_vola_forecast(1:foptions.S), 'k--','LineWidth',2);
plot(dates(cutoff-49:cutoff), tau_forecast_annual(end-49-foptions.S:end-foptions.S),'r-', 'LineWidth',2);
plot(dates(cutoff+1:cutoff+foptions.S), tau_forecast_annual(end-foptions.S+1:end),'r--','LineWidth',2);
yline(annual_unconditional_vola, 'b--'); 
xlabel('Date', 'Interpreter', 'latex'); ylabel('conditional/long-term volatility', 'Interpreter', 'latex'); hold off 
xlim([min(dates(cutoff-49:cutoff+foptions.S)) max(dates(cutoff-49:cutoff+foptions.S))]);
xline(cutoff_date+1, '--', 'LineWidth', 1);
xl.LabelVerticalAlignment = 'middle';
xl.LabelHorizontalAlignment = 'center';
legend({'$\sigma_t$ conditional volatility (in-sample)', ...
        '$\sigma_t$ conditional volatility (out-of-sample forecast)', ... 
        '$\tau_t$ long-term volatility (in-sample)', ...
        '$\tau_t$ long-term volatility (out-of-sample forecast)', ...
        'unconditional volatility'
        }, ...
        'Interpreter', 'latex', 'Location', 'best');
set(gcf, 'Color', 'w');
set(gca,'TickLabelInterpreter','latex', 'FontSize', 14); 
saveas(gcf,'figures/ForecastIllustration.png'); 
end 