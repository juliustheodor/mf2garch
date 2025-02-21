function [ll]=likelihood_mf2_garch(param, y, m)

% Input: 
%   parameters: (7x1) coefficient vector for the MF2-GARCH, estimation results for
%       (mu, alpha, gamma, beta, lambda_0, lambda_1, lambda_2)' 
%   y: a (Tx1) vector with the return time series 
%   m: a scalar (used for the rolling window measure of the local bias of 
%      the short-term component’s conditional variance over the previous m days) 

% Output: 
    % ll: Log-likelihood-function (scalar) 

    [e, h, tau, V_m ] = mf2_garch_core(param, y, m);

    lls = -0.5*(log(2*pi) + log(h.*tau) + e.^2);
    
    % Use these to comput the LL
    ll = -sum(lls);
    
end