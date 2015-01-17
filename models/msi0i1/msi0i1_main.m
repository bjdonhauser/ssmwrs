function results = msi0i1_main()
% msi0i1_main  --  The main method for the msi0i1 model. As presented in Kim 
%               and Nelson (1999, pg. 111). Returns parameter estimates, 
%               standard errors, state inferences, state variable 
%               inferences, likelihood value, and other output.
%
%
%****f* SSMWRS/models/msi0i1/msi0i1_main
%
% NAME
%   msi0i1_main --  The main method for the msi0i1 model. Returns parameter 
%                   estimates, standard errors, state inferences, state 
%                   variable inferences, likelihood value, and other 
%                   output.
%
% SYNOPSIS
%   [results] = msi0i1_main()
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
%   The general form for a state-space model with regime switching 
%   (following the notation of Kim and Nelson (1999, pg. 98)) is the 
%   following:
%
%   Measurement equation:
%   y_t = H_S_t*b_t + A_S_t*z_t + e_t           e_t ~ N(0, R_S_t)
%
%   Transition equation:
%   b_t = mu_S_t + F_S_t*b_(t-1) + G_S_t*v_t    v_t ~ N(0, Q_S_t)
%
%   Transitional probability matrix:
%   p = [p11, p21, ... , pM1; p12, ..., pM2; ... ; p1M, ... pMM]; where pjk
%   = Pr{S_t=k|S_(t-1)=j}
%
%   The msi0i1 model is intended to be an alternative to the Stock and
%   Watson (2006) UC-SV model of inflation. In the msi0i1 model there is
%   Markov switching between inflation having an I(0) and I(1) component.
%   The model looks as follows:
%
%   pi(t) = (1)*tau(t) + eta            eta ~ N(0,sig2eta)
%   tau(t) = (1)*tau(t-1) + epsilon     epsilon ~ N(0,sig2eps_j)
%
%   sigma_epsilon_j^2 = sigma_epsilon_1^2 when S_t = 2 and 0 when S_t = 1
%
%   Pr{S_t=2|S_(t-1)=2}=p22, Pr{S_t=1|S_(t-1)=1}=p11
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   msi0i1_constraints, msi0i1_ssf, ssm_optimizer, ssm_objective,
%   kim_filter.
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-12
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
parms_0 = [1,0.5,5,5]';  
model_ssf = @msi0i1_ssf;
model_constraints = @msi0i1_constraints;
opt_engine = @fminunc;
logL_start = 25;    % 1959:1. Forecasting will also begin for 1959:1. 
                    % Follows Stock and Watson
opt = optimset('TolX', 0.0001, 'Display', 'off', 'MaxIter', 200, ...
               'MaxFunEvals', 200, 'LargeScale', 'off');
filter = @kim_filter;

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
okh_tsplot(extract_results_ts(results.S_t_j,2,logL_start),info);
title('Probability of being in state 2 (the I(1) state)');

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



%=========================================================================%
% Results
%=========================================================================%

%===============================%
% #1
%===============================%
% parms_0 = [1,0.5,5.0,5.0]'; % corresponds to [1,0.25,0.9615,0.9615]' 
% parms_est: [0.7005 1.1405 0.9193 0.9324]
% parms_se: [0.13416;0+i*0.84198;0+i*0.045523;0+i*0.20718]
% logL: -302.1414
% MSFE: 1.43    (from:
%               foo = extract_results_ts(results.h_l,1,25);
%               sum(foo.^2)/193
% COMMENTS: (1) transitory variance smaller. (2) Filtered estimate of the
% state pretty erratic. Although, we atleast get values that we expect.

%===============================%
% #2
%===============================%
% parms_0 = [1,0.5,0.5,0.5]'; % corresponds to [1,0.25,0.20,0.20]' 
% parms_est: [0.5507 2.3054 0.8257 0.5999]
% parms_se: [0.11172;1.5296;0.070564;0.38959]
% logL: -300.6408
% MSFE: 1.45915
% COMMENTS: (1) A similar fit to the first. Maybe more realistic

%===============================%
% #3
%===============================%
% parms_0 = [1,0.5,1,1]'; % corresponds to [1,0.25,0.50,0.50]' 
% parms_est: [0.5450 2.1427 0.8336 0.5592]
% parms_se: [0.10785;1.1209;0.07312;0.3488]
% logL: -300.5744
% MSFE: 1.4580

%===============================%
% #4
%===============================%
% parms_0 = [1,0.5,2,2]'; % corresponds to [1,0.25,0.80,0.80]' 
% parms_est: [0.5348 1.6332 0.8621 0.7639]
% parms_se: [0.11389;0.69412;0.070501;0.18466]
% logL: -300.1634
% MSFE: 1.4527

%===============================%
% #5: NAIVE FORECAST
%===============================%
% MSFE: 1.6986
