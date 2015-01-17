function parms_constrained = lam_constraints(parms)
% lam_constraints -- Transforms the unconstrained parameters of the lam 
%                    model to constrained parameters. This function is
%                    intended to be called by ms_ssm_objective. The purpose
%                    of this function is so that unconstrained optimization
%                    techniques may be used (which sounds counterintuitive).
%
%****f* SSMWRS/lam_constraints
%
% NAME
%   lam_constraints --  Transforms the unconstrained parameters of the lam 
%                       model to constrained parameters. This function is
%                       intended to be called by ms_ssm_objective. 
%
% SYNOPSIS
%   [SSF] = lam_constraints(parms)
%
% INPUTS
%   * parms             -- vector of length 7 corresponding to [p, q,
%                          delta_0, delta_1, sigma, phi_1, phi_2]'.
%
% OUTPUTS
%   * parms_constrained -- the constrained parameters
%
% SIDE EFFECTS
%
% DESCRIPTION
%   The constraints are intended to keep p and q positive and less than 1,
%   sigma positive, and phi_1 and phi_2 satisfying stationary values. 
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   lam_main, lam_ssf, kim_filter, ms_ssm_optimizer, ms_ssm_objective
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-10
%
%***

% Initialize
parms_constrained = zeros(7,1);

% Constrains so that 0 < p, q < 1
% parms_constrained(1:2) = parms(1:2).^2./(1+parms(1:2).^2);

% Using the exact kim_je.opt method instead.
parms_constrained(1:2) = exp(parms(1:2))./(1+exp(parms(1:2)));

% No constraints on delta_0 or delta_1
parms_constrained(3:4) = parms(3:4);

% Constrains so that sigma > 0
% NOTE: don't really have to constrain this because it's squared in the
% model
% parms_constrained(5) = parms(5)^2;
parms_constrained(5) = parms(5);

% Constrains so that phi_1 and phi_2 jointly satisfy stationarity
% properties (see Hamilton (1994))
temp = parms(6:7)./(1+abs(parms(6:7)));
parms_constrained(6) = sum(temp);
parms_constrained(7) = -1*temp(1)*temp(2);


