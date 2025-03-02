function [coeff, qmle_se, p_value_qmle, Z, h, tau, sigma_annual, tau_annual, annual_unconditional_vola, foptions] = mf2_garch_estimation(y, foptions)

% MF2 GARCH-rw-m Toolbox for Matlab by Christian Conrad & Julius Schoelkopf
% Version 0.1.0

% INPUTS: 
%   y: a (Tx1) vector with the return time series 
%   foptions: 
%       f.options.m: a scalar (used for the rolling window measure of the local bias of 
%       the short-term component’s conditional variance over the previous m days) 
%       foptions.choice = 'BIC' or 'fix'
%           - 'BIC': choose the optimal m as the one that minimizes the BIC   
%           - 'fix': you need to specify the m you want to use using
%               f.options.m. 

% OUTPUTS: 
%   coeff: (7x1) coefficient vector for the MF2-GARCH, estimation results for
%       (mu, alpha, gamma, beta, lambda_0, lambda_1, lambda_2)' 
%   qmle_se: (7x1) vector with Bollerslev–Wooldridge robust standard errors 
%   p_value_qmle: (7x1) vector with p-values based on Bollerslev–Wooldridge 
%     robust standard errors 
%   BIC: (Normalized) Bayesian information criterion
%   Z: (length(y) - 2*252 x 1 ) vector with (e= (y-mu) ./ sqrt(h.*tau)) 
%   h: (length(y) - 2*252 x 1 ) vector of fitted values for short-term component 
%   tau: (length(y) - 2*252 x 1 ) vector fitted values for long-term component 
%   sigma_annual: (length(y) - 2*252 x 1 ) vector fitted values for annualized cond. volatility 
%   tau_annual: (length(y) - 2*252 x 1 ) vector fitted values for annualized long-term component 
%   annual_unconditional_vola: scalar of annualized unconditional volatility 
%   foptions: extended with the optimal m according to the BIC if BIC
%       selection was specified. 

% Number of observations
t = length(y); 

% parameters: mu, alpha, gamma, beta, lambda_0, lambda_1, lambda_2
% start values
param_init = [0.02; 0.007; 0.14 ; 0.85; mean(y.^2)*(1-0.07-0.91); 0.07; 0.91 ];

% constraints:
% alpha >=0
% alpha + gamma/2 + beta <=1
% lambda_1 >= 0
% lambda_1 + lambda_2 <=1
A = [  ...
    0.0 -1.0  0.0  0.0  0.0  0.0  0.0; ...
    0.0 1.0  0.5  1.0  0.0  0.0  0.0; ...
    0.0 0.0  0.0  0.0  0.0  -1.0  0.0; ...
    0.0 0.0  0.0  0.0  0.0  1.0  1.0...
];
b = [  0.0; 1.0; 0.0; 1.0 ];

Aeq = [];
beq = [];
nlc = [];

% bounds:
LB = [-1; 0.0 ; -0.5; 0.0; 0.000001; 0.0; 0.0 ];
UB = [ 1; 1.0 ; 0.5; 1.0; 10.0; 1.0; 1.0 ];

options = optimoptions(@fmincon,'Algorithm','sqp'); 


if isfield(foptions, 'choice') && strcmp(foptions.choice, 'BIC')

    BIC_vec=zeros(130,1);  

for cm = 20:150
    [~, ll, ~, ~, ~, ~, ~] = fmincon('likelihood_mf2_garch', param_init, A, b, Aeq, beq, LB, UB, nlc, options, y, cm); 
    
    % Information Criteria 
    [~,~,ic] = aicbic(-ll,7,t, Normalize=true); 
    BIC_vec(cm-19) = ic.bic; 

end 

minimum =  min(BIC_vec);
m = find(BIC_vec ==minimum)+19; 

fprintf('Optimal m =  %d', m);

foptions.m = m; 

% Plot of BIC
figure
plot(20:150, BIC_vec, 'k-'); hold on;
xline(m, 'r-', 'LineWidth', 1); hold off 
xlabel('$m$', 'Interpreter', 'latex'); ylabel('BIC', 'Interpreter', 'latex'); 
set(gcf, 'Color', 'w');
set(gca,'TickLabelInterpreter','latex', 'FontSize', 14); 
saveas(gcf,'figures/BIC.png'); 

else
    m = foptions.m; 
end


%% Estimation of the MF2-GARCH-rw-m (see section A.1.1 in Conrad & Engle (2025) for details) 
[ coeff, ll, exitFlag, ~, ~, grad, hessian ] = fmincon('likelihood_mf2_garch', param_init, A, b, Aeq, beq, LB, UB, nlc, options, y, m); 
[ Z, h, tau, V_m ] = mf2_garch_core(coeff, y, m);
[qmle_se, p_value_qmle] = mf2_garch_inference(coeff, y, Z, h, tau, m); 

% Annualize Data (required for forecasting plots) 
sigma_annual = sqrt(252.*tau.*h);

% Annualize conditional volatility   
tau_annual = sqrt(252*tau);

%% Display Results in Table 
% Example 7x1 vectors
labels = {'mu'; 'alpha'; 'gamma'; 'beta'; 'lambda_0'; 'lambda_1'; 'lambda_2'}; 

% Number of estimated parameters
k = length(coeff); 

% Information Criteria 
[~,~,ic] = aicbic(-ll,7,t, Normalize=true); 

% Kappa
kappa = mean((Z.^2 - 1).^2) + 1; 

% Check covariance-stationarity 

mu = coeff(1);
alpha = coeff(2);
gamma = coeff(3);
beta = coeff(4);
    
lambda_0 = coeff(5);
lambda_1 = coeff(6);
lambda_2 = coeff(7);

[delta_m] = Delta_m(alpha, gamma, beta, lambda_0, lambda_1, lambda_2, kappa, m);

numerator = lambda_0 + lambda_0/(1-lambda_1-lambda_2) * (1-(alpha + gamma/2)-beta) * (lambda_1 + lambda_2) + delta_m;

[sum_alpha_beta_deno] = alpha_beta_deno(alpha, gamma, beta, m);
    
Gamma_m = (lambda_1*1/m*((alpha + gamma/2)*kappa+beta) + lambda_2*((alpha + gamma/2)+beta)) - lambda_1*1/m*((alpha + gamma/2)*kappa+beta) * sum_alpha_beta_deno;

denominator = 1 - (lambda_1*1/m*((alpha + gamma/2)*kappa+beta) + lambda_2*((alpha + gamma/2)+beta)) - lambda_1*1/m*((alpha + gamma/2)*kappa+beta) * sum_alpha_beta_deno;
	
var_m = numerator / denominator;
     
annual_unconditional_vola = (252 * var_m)^0.5; 

% Define significance levels
significance = strings(size(p_value_qmle')); % Preallocate an empty string array

% Assign stars based on p-value thresholds
significance(p_value_qmle' < 0.01) = "***";
significance(p_value_qmle' >= 0.01 & p_value_qmle' < 0.05) = "**";
significance(p_value_qmle' >= 0.05 & p_value_qmle' < 0.10) = "*"; 

% Create a table
T = table(labels, coeff, qmle_se', p_value_qmle', significance,...
    'VariableNames', {'Parameter', 'Coefficient', 'Standard Error', 'p-value', 'Significance'});


% Display the table
fprintf('\n===================== Estimation results MF2-GARCH-rw-m =====================\n\n');
if isfield(foptions, 'choice') && strcmp(foptions.choice, 'BIC')
     fprintf('The optimal m was selected as the one that minimizes the BIC: m = %d\n', m); 
else 
     fprintf('The optimal m was specified by the user: m = %d\n', m); 
end 
fprintf('Log-Likelihood Function = %.3f, BIC = %.3f\n', -ll, ic.bic);
fprintf('Estimated fourth moment of the innovations: kappa = %.3f\n', kappa);

disp(T)
fprintf('Output reports Bollerslev-Wooldridge robust standard errors (see Conrad and \nEngle (2025), equation (27)).\n'); 
if (Gamma_m < 1)
    fprintf('Covariance stationarity condition satisfied (see Conrad and Engle (2025), \nequation (7)): Gamma_m = %.3f\n', Gamma_m);    
else
    fprintf('Covariance stationarity condition violated (see Conrad and Engle (2025),\nequation (7)): Gamma_m = %.3f\n',  Gamma_m);    
end
fprintf('Annualized unconditional volatility = %.3f\n', annual_unconditional_vola);
fprintf('==============================================================================\n');

foptions.BIC = ic.bic; 
end 