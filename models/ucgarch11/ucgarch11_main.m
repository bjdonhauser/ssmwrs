function results = ucgarch11_main()
% ucgarch11_main  --  The main method for the ucgarch11 model. As presented in Kim 
%               and Nelson (1999, pg. 111). Returns parameter estimates, 
%               standard errors, state inferences, state variable 
%               inferences, likelihood value, and other output.
%
%
%****f* SSMWRS/models/ucgarch11/ucgarch11_main
%
% NAME
%   ucgarch11_main --  The main method for the ucgarch11 model. Returns parameter 
%                   estimates, standard errors, state inferences, state 
%                   variable inferences, likelihood value, and other 
%                   output.
%
% SYNOPSIS
%   [results] = ucgarch11_main()
%
% INPUTS
%
% OUTPUTS
%   * results --    parameter estimates, standard errors, state 
%                   inferences, state variable inferences, likelihood 
%                   value, and other output. (structure)
%
% SIDE EFFECTS
%
% DESCRIPTION
%   The general form for a state-space model with GARCH errors is given in
%   Kim and Nelson (1999, pg. 140) and the augmented form pgs. 142-143.
%   This is intended to mimic the UC-SV model of Stock and Watson. Where we
%   have:
%
%   Measurement Equation:
%   pi_t = tau_t + eta_t, eta_t ~ GARCH(1,1)
%
%   Transition Equation:
%   tau_t + tau_(t-1) + epsilon_t, epsilon_t ~ GARCH(1,1)
%
%   eta_t and epsilon_t should be more or less independently distributed.
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   ucgarch11_constraints, ucgarch11_ssf, kalman_filter_garch.
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-13
%
%***
%=========================================================================%
% Load and process data
%=========================================================================%
PCE = xlsread('Z:/data/PCE/PCECTPI');               %1947:1-2007:1
inflation = calc_growthrate(PCE, 'quarterly');      %1947:2-2007:1        
inflation = inflation(24:length(inflation));        %1953:1-2007:1 (follows
                                                    %Stock and Watson
N = length(inflation);

%=========================================================================%
% Set parameters 
%=========================================================================%
parms_0 = [-1, 0.5, 0.5, -1, 0.5, 0.5]';            % [alpha_0, alpha_1, 
                                                    % alpha_2, gamma_0, 
                                                    % gamma_1, gamma_2]
model_ssf = @ucgarch11_ssf;
model_constraints = @ucgarch11_constraints;
opt_engine = @fminsearch;
logL_start = 25;    % 1959:1. Forecasting will also begin for 1959:1. 
                    % Follows Stock and Watson
opt = optimset('TolX', 0.0001, 'Display', 'on', 'MaxIter', 5000, ...
               'MaxFunEvals', 5000, 'LargeScale', 'on');
filter = @kalman_filter_garch;

%=========================================================================%
% Call the optimizer
%=========================================================================%
results = ssm_optimizer(parms_0, inflation, model_ssf, model_constraints, ...
                           opt_engine, logL_start, opt, filter);

%=========================================================================%
% Plotting
%=========================================================================%

figure;
info.yy   = 1959;
info.pr   = 1;
info.freq = 'Q';
okh_tsplot(extract_results_ts(results.b_t,1,logL_start),info);
hold on;
okh_tsplot(inflation(logL_start:N),info,'Color','r');
title('Filtered estimate of the trend (inflation series in red)');
hold off;

figure;
info.yy   = 1959;
info.pr   = 1;
info.freq = 'Q';
okh_tsplot(extract_results_ts(results.y_l,1,logL_start),info);
hold on;
okh_tsplot(inflation(logL_start:N),info,'Color','r');
title('Filtered 1-step ahead forecast of inflation (inflation series in red)');
hold off;

figure;
info.yy   = 1959;
info.pr   = 1;
info.freq = 'Q';
okh_tsplot(extract_results_ts(results.eta_l,1,logL_start),info);
title('1-step ahead Prediction Errors');

figure;
info.yy   = 1959;
info.pr   = 1;
info.freq = 'Q';
okh_tsplot(extract_results_ts(results.f_l,1,logL_start),info);
title('Conditional variance of 1-step ahead Prediction Errors');

figure;
info.yy   = 1959;
info.pr   = 1;
info.freq = 'Q';
okh_tsplot(extract_results_ts(results.h_me_t,1,logL_start),info);
title('Conditional Variance of the Transitory Distrubance');

figure;
info.yy   = 1959;
info.pr   = 1;
info.freq = 'Q';
okh_tsplot(extract_results_ts(results.h_te_t,1,logL_start),info);
title('Conditional Variance of the Permanent Distrubance');



%=========================================================================%
% Results
%=========================================================================%
