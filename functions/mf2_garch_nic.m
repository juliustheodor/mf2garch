function [r, NIC ] = mf2_garch_nic(Z, h, tau, foptions, coeff)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

% INPUTS: 
%   Z: (length(y) - 2*252 x 1 ) vector with (e= (y-mu) ./ sqrt(h.*tau)) 
%   h: (length(y) - 2*252 x 1 ) vector of fitted values for short-term component 
%   tau: (length(y) - 2*252 x 1 ) vector fitted values for long-term component 
%   foptions: 
%       f.options.m: a scalar (used for the rolling window measure of the local bias of 
%       the short-term component’s conditional variance over the previous m days) 
%   coeff: (7x1) coefficient vector for the MF2-GARCH, estimation results for
%       (mu, alpha, gamma, beta, lambda_0, lambda_1, lambda_2)' 

% OUTPUTS: 
%   NIC: news-impact curve 
%   r: return grid 

    % parameters 
    mu = coeff(1);
    alpha = coeff(2);
    gamma = coeff(3);
    beta = coeff(4);
    
    lambda_0 = coeff(5);
    lambda_1 = coeff(6);
    lambda_2 = coeff(7);
    
    m = foptions.m; 
    
    steps = 100;
    last = steps + 1;
    
    r = -5.1 .* ones(last,1);
    h_t_1 = zeros(last,1);
    tau_t_1 = zeros(last,1);
    
    NIC = zeros(last,1);
    
    V = Z.^2 .* tau;
   
    const = 0.1;
    
    for t = 1:last
   
    r(t) = r(t) + const;
       
    if ((r(t) - mu) < 0) 
        
    h_t_1(t) = (1-alpha-gamma/2-beta) + (alpha + gamma).*(r(t) - mu).^2 ./ tau(length(Z)) + beta .* h(length(Z));
    
    else
        
    h_t_1(t) = (1-alpha-gamma/2-beta) + alpha .* (r(t) - mu).^2 ./ tau(length(Z)) + beta .* h(length(Z));
    
    end
       
    tau_t_1(t) = lambda_0 + lambda_1 ./m * (r(t) - mu).^2 ./ h(length(Z)) + lambda_1 ./m * sum(V(length(Z)-(m-1):length(Z)-1)) + lambda_2 * tau(length(Z));
    
    NIC(t) = (252 .* h_t_1(t) .* tau_t_1(t)).^0.5;
    
    const = const + 0.1;
    
    end

    % Figure 
    figure
    plot(r, NIC, 'k','LineWidth',2);
    title('News Impact Curve for MF2-GARCH-rw-m', 'Interpreter', 'latex', 'FontSize', 12); 
    legend(sprintf('m = %d', m), 'Interpreter','latex');
    set(gcf, 'Color', 'w');
    xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5]);  
    xline(0, 'k', 'HandleVisibility', 'off');
    set(gca,'TickLabelInterpreter','latex', 'FontSize', 14); 
    ylim([0 inf])
    saveas(gcf,'figures/NIC.png'); 

end 
   
    