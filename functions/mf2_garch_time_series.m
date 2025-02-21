function  [] = mf2_garch_time_series(dates, sigma_annual, tau_annual)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

% INPUTS: 
%   dates: date vector (could also be replaced with a 1:T vector) 
%   sigma_annual: (length(y) - 2*252 x 1 ) vector fitted values for annualized cond. volatility 
%   tau_annual: (length(y) - 2*252 x 1 ) vector fitted values for annualized long-term component 

% OUTPUT: 
% Figure saved as 'TimeSeries.png' in the figures folder.

    figure
    plot(dates(505:end), sigma_annual,'k-','LineWidth',1); hold on; 
    plot(dates(505:end), tau_annual,'r-','LineWidth',2); 
    recessionplot; 
    legend({'$\sigma_t$ conditional volatility', ...
            '$\tau_t$ long-term volatility', 
            }, ...
            'Interpreter', 'latex', 'Location', 'best');
    set(gcf, 'Color', 'w');
    set(gca,'TickLabelInterpreter','latex', 'FontSize', 14); 
    saveas(gcf,'figures/TimeSeries.png'); 
end 