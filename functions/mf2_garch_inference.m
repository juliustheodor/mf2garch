% Calculation of Bollerslev–Wooldridge robust standard errors 

function [qmle_se, p_value_qmle] = mf2_garch_inference(parameters, y, e, h, tau, m)

% INPUTS: 
%   parameters: (7x1) coefficient vector for the MF2-GARCH, estimation results for
%       (mu, alpha, gamma, beta, lambda_0, lambda_1, lambda_2)' 
%   qmle_se: (7x1) vector with Bollerslev–Wooldridge robust standard errors 
%   p_value_qmle: (7x1) vector with p-values based on Bollerslev–Wooldridge 
%     robust standard errors 
%   y: a (Tx1) vector with the return time series 
%   e: (length(y) - 2*252 x 1 ) vector with (e= (y-mu) ./ sqrt(h.*tau)) 
%   h: (length(y) - 2*252 x 1 ) vector of fitted values for short-term component 
%   tau: (length(y) - 2*252 x 1 ) vector fitted values for long-term component 
%   m: a scalar (used for the rolling window measure of the local bias of 
%      the short-term component’s conditional variance over the previous m days) 

% OUTPUTS: 
%   qmle_se: (7x1) vector with Bollerslev–Wooldridge robust standard errors 
%   p_value_qmle: (7x1) vector with corresponding p-values based on 
%       Bollerslev–Wooldridge robust standard errors 

% This code uses a function for 2-sided finite difference Hessian from the 
% Oxford MFE Toolbox from Kevin Sheppard, Available from https://github.com/bashtage/mfe-toolbox

    mu = parameters(1);
    alpha = parameters(2);
    gamma = parameters(3);
    beta = parameters(4);
    
    lambda_0 = parameters(5);
    lambda_1 = parameters(6);
    lambda_2 = parameters(7);

    T = length(tau); 
    k=length(parameters);

    %log-likelihood contributions at estimates 
    lls = -0.5*(log(2*pi) + log(h.*tau) + e.^2);
    
    % Two-sided finite Differences for each parameter 
    hhh=max(abs(parameters'*eps^(1/3)),1e-8);
    hhh = diag(hhh); 
    scores=zeros(T,k);

    for j=1:k %for each parameter

        parameters_h_p = parameters + hhh(:,j);
        
        [ e_h, h_h, tau_h, ~] = mf2_garch_core(parameters_h_p, y, m);

        lls_p = -0.5*(log(2*pi) + log(h_h.*tau_h) + e_h.^2);
        
        parameters_h_m = parameters - hhh(:,j); 

        [ e_m, h_m, tau_m, ~] = mf2_garch_core(parameters_h_m, y, m);

        lls_m = -0.5*(log(2*pi) + log(h_m.*tau_m) + e_m.^2);

        %using log-likelihood contributions 
        scores(:,j) = (lls_p - lls_m)./(2.*hhh(j,j)) ; 
        
    end 
    
    S = (1/T).*scores'*scores; 

    % Hessian 
    H = hessian_2sided('likelihood_mf2_garch',parameters, y,m); 
    A = H/T; 
    mhess = A^(-1); 
    
    % Standard Errors 
    qmle_se = sqrt(diag((mhess*(S)*mhess)/T))'; 

    % p-value 
    p_value_qmle = 2*(1-  normcdf(abs(parameters'./qmle_se),0,1)) ; 
end 
