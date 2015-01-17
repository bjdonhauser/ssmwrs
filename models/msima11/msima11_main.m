function results = msima11_main()
% msima11_main  --  The main method for the msima11 model. As presented in Kim 
%               and Nelson (1999, pg. 111). Returns parameter estimates, 
%               standard errors, state inferences, state variable 
%               inferences, likelihood value, and other output.
%
%
%****f* SSMWRS/models/msima11/msima11_main
%
% NAME
%   msima11_main --  The main method for the msima11 model. Returns parameter 
%                   estimates, standard errors, state inferences, state 
%                   variable inferences, likelihood value, and other 
%                   output.
%
% SYNOPSIS
%   [results] = msima11_main()
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
%   The msima11 model is intended to be an alternative to the Stock and
%   Watson (2006) UC-SV model of inflation. In the msima11 model there is
%   Markov switching between inflation having an I(0) and I(1) component.
%   by virtue of how theta switches from theta = 1 to other values:
%
%   diff( pi(t) )= epsilon_t - theta_S_t*epsilon_t     epsilon_t ~ N(0,sig2)
%   tau(t) = (1)*tau(t-1) + epsilon     epsilon ~ N(0,sig2eps_j)
%
%   theta_S_t   = 1         for S_t = 1
%               = theta_2   for S_t = 2
%
%   Pr{S_t=2|S_(t-1)=2}=p22, Pr{S_t=1|S_(t-1)=1}=p11
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   msima11_constraints, msima11_ssf, ssm_optimizer, ssm_objective,
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
dInflation = diff(inflation);                       %1947:3-2007:1
inflation = inflation(23:length(inflation));        %1953:1-2007:1 (follows
                                                    %Stock and Watson

%=========================================================================%
% Set parameters 
%=========================================================================%
parms_0 = [1,0.5,2,2]';  
model_ssf = @msima11_ssf;
model_constraints = @msima11_constraints;
opt_engine = @fminunc;
logL_start = 24;    % 1959:1. Forecasting will also begin for 1959:1. 
                    % Follows Stock and Watson
opt = optimset('TolX', 0.0001, 'Display', 'off', 'MaxIter', 50, ...
               'MaxFunEvals', 50, 'LargeScale', 'off');
filter = @kim_filter;

%=========================================================================%
% Call the optimizer
%=========================================================================%
results = ssm_optimizer(parms_0, dInflation, model_ssf, model_constraints, ...
                           opt_engine, logL_start, opt, filter);


%=========================================================================%
% Results
%=========================================================================%

%===============================%
% #1
%===============================%
% parms_0 = [1,0.5,2,2]'; % corresponds to [1,0.25,0.8,0.8]' 
% parms_est: [0.0510 1.0741 0.8576 0.7006]
% parms_se: [0.087663;0.1129;0.10587;0.18364]
% logL: -351.77
% MSFE: 1.5066    (from:
%               foo = extract_results_ts(results.eta_l,1,25);
%               sum(foo.^2)/193

