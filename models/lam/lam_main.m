function results = lam_main()
% lam_main  --  The main method for the lam model. As presented in Kim 
%               and Nelson (1999, pg. 111). Returns parameter estimates, 
%               standard errors, state inferences, state variable 
%               inferences, likelihood value, and other output.
%
%
%****f* SSMWRS/models/lam/lam_main
%
% NAME
%   lam_main    --  The main method for the lam model. As presented in Kim 
%                   and Nelson (1999, pg. 111). Returns parameter 
%                   estimates, standard errors, state inferences, state 
%                   variable inferences, likelihood value, and other 
%                   output.
%
% SYNOPSIS
%   [results] = lam_main()
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
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   lam_constraints, lam_ssf, ssm_optimizer, ssm_objective,
%   kim_filter.
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-12
%
%***

% NOTE: When the starting values and parameter estimates of the book are
% used, the results are not replicated. However, when the starting values
% and parameter estimates of the program kim_je.opt are used, we get
% perfect replication. There must be more happening behind the scenes of
% table 5.1 in Kim and Nelson.

%=========================================================================%
% Load and process data
%=========================================================================%
GNP = xlsread('Z:/data/GDP/kim_je.xls');
gGNP = calc_growthrate(GNP);                                %1947:2-1986:4   

N = length(gGNP);
gGNP = gGNP(1:(N-8));                                       %1947:2-1984:4
N = length(gGNP);

%=========================================================================%
% Set parameters 
%=========================================================================%
parms_0 = [3.5, 0.0, -0.4, 1.3, 0.7, 0.5, 1]'; 
model_ssf = @lam_ssf;
model_constraints = @lam_constraints;
opt_engine = @fminunc;
logL_start = 23;
opt = optimset('TolX', 0.0001, 'Display', 'off', 'MaxIter', 5000, ...
               'MaxFunEvals', 5000, 'LargeScale', 'off');
filter = @kim_filter;

% set parameters as the final parameters of kim_je.opt
% NOTE!!! While my delta_1 = their delta_0, my delta_2 = their delta_0 +
% delta_1!!! This difference shows up in parms(4).
%parms = [0.9503,0.4428,-1.2917,0.9458,0.8014,1.2608,-0.3534]';

%=========================================================================%
% Call the optimizer
%=========================================================================%
results = ssm_optimizer(parms_0, gGNP, model_ssf, model_constraints, ...
                           opt_engine, logL_start, opt, filter);

