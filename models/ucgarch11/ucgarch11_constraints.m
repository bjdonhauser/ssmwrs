function parms_constrained = ucgarch11_constraints(parms)
% ucgarch11_constraints --  Transforms the unconstrained parameters of 
%                           the ucgarch11 model to constrained parameters. 
%                           The purpose of doing this is so that and 
%                           unconstrained optimization can be performed in
%                           the new, transformed variables
%
%****f* SSMWRS/models/ucgarch11/ucgarch11_constraints
%
% NAME
%   ucgarch11_constraints --	Transforms the unconstrained parameters of 
%                               the ucgarch11 model to constrained 
%                               parameters. The purpose of doing this is 
%                               so that and unconstrained optimization can 
%                               be performed in the new, transformed 
%                               variables
%
% SYNOPSIS
%   [SSF] = ucgarch11_constraints(parms)
%
% INPUTS
%   * parms             -- vector of length 6 corresponding to the values
%                          [alpha_0, alpha_1, alpha_2, gamma_0, gamma_1, 
%                          gamma_2]
%
% OUTPUTS
%   * parms_constrained -- the constrained parameters
%
% SIDE EFFECTS
%
% DESCRIPTION
%   Must constrain GARCH parameters so that the variance of the distrubance
%   term follows a stationary process. For GARCH(1,1) this means that
%   alpha_0 > 1 and 0 < alpha_1 + alpha_2 < 1.
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   ucgarch11_main, ucgarch11_ssf, kalman_filter_garch
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-13
%
%***

%=========================================================================%
% Initialize
%=========================================================================%
parms_constrained = zeros(6,1);

%=========================================================================%
% Constraining Parameters for the GARCH distrubance term in the Measurement
% equation
%=========================================================================%

% alpha_0 > 1
parms_constrained(1) = exp(parms(1)); 

% 0 < alpha_1 + alpha_2 < 1;
parms_constrained(2) = exp(parms(2))/(1+exp(parms(2))+exp(parms(3)));
parms_constrained(3) = exp(parms(3))/(1+exp(parms(2))+exp(parms(3)));

%=========================================================================%
% Constraining Parameters for the GARCH distrubance term in the Transition
% equation
%=========================================================================%

% gamma_0 > 1
parms_constrained(4) = exp(parms(4)); 

% 0 < gamma_1 + gamma_2 < 1;
parms_constrained(5) = exp(parms(5))/(1+exp(parms(5))+exp(parms(6)));
parms_constrained(6) = exp(parms(6))/(1+exp(parms(5))+exp(parms(6)));


