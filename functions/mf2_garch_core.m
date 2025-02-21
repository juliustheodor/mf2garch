function [ e, h, tau, V_m ] = mf2_garch_core(parameters, y, m)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

% INPUTS: 
%   parameters: (7x1) coefficient vector for the MF2-GARCH, estimation results for
%       (mu, alpha, gamma, beta, lambda_0, lambda_1, lambda_2)' 
%   y: a (Tx1) vector with the return time series 
%   m: a scalar (used for the rolling window measure of the local bias of 
%      the short-term component’s conditional variance over the previous m days) 

% OUTPUTS: 
%   e: (length(y) - 2*252 x 1 ) vector with (e= (y-mu) ./ sqrt(h.*tau)) 
%   h: (length(y) - 2*252 x 1 ) vector of fitted values for short-term component 
%   tau: (length(y) - 2*252 x 1 ) vector fitted values for long-term component 
%   V_m: rolling window measure of the local bias of the short-term
%       component’s conditional variance over the previous m days.

    mu = parameters(1);
    alpha = parameters(2);
    gamma = parameters(3);
    beta = parameters(4);
    
    lambda_0 = parameters(5);
    lambda_1 = parameters(6);
    lambda_2 = parameters(7);

    h = ones(size(y));
    
    tau = ones(size(y))*mean(y.^2);
    V = zeros(size(y));
    V_m = zeros(size(y));
    
    
    
    for t = 2:m
        
       if ((y(t-1) - mu) < 0) 
       h(t) = (1-alpha-gamma/2-beta) + (alpha + gamma).*(y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1);
       else
       h(t) = (1-alpha-gamma/2-beta) + alpha .* (y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1); 
       end
       
    end
    
    
    for t = (m+1):length(y)
        
       if ((y(t-1) - mu) < 0) 
       h(t) = (1-alpha-gamma/2-beta) + (alpha + gamma).*(y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1);
       else
       h(t) = (1-alpha-gamma/2-beta) + alpha .* (y(t-1) - mu).^2 ./ tau(t-1) + beta .* h(t-1); 
       end
       
       V(t) = (y(t) - mu).^2 ./ h(t);
       V_m(t) = sum(V(t-(m-1):t))./m;
       
       tau(t) = lambda_0 + lambda_1 * V_m(t-1) + lambda_2 * tau(t-1);
       
    end    
    
    e = (y-mu) ./ sqrt(h.*tau);
    
    sample = (2*252+1):length(y);
    h=h(sample);
    tau=tau(sample);
    e=e(sample);
    
    end

