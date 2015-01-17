function parms_constrained = msima11_constraints(parms)
% msima11_constraints -- Transforms the unconstrained
% parameters of the msima11 model to constrained parameters. The purpose of
% doing this is so that and unconstrained optimization can be performed in
% the new, transformed variables
%
%****f* SSMWRS/models/msima11/msima11_constraints
%
% NAME
%   msima11_constraints --  Transforms the unconstrained
% parameters of the msima11 model to constrained parameters. The purpose of
% doing this is so that and unconstrained optimization can be performed in
% the new, transformed variables  
%
% SYNOPSIS
%   [SSF] = msima11_constraints(parms)
%
% INPUTS
%   * parms             -- vector of length 4 corresponding to the values
%                          [theta_2, sig2, p11, and p22]
%
% OUTPUTS
%   * parms_constrained -- the constrained parameters
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
%   ssm_optimizer, ssm_objective, msima11_main, msima11_ssf
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-07
%
%***

% NOTE: the constraining method avoids using exponentials (as in
% Clark_rest.m) which are more likely to lead to problems for the
% optimizer. 
parms_constrained(1) = parms(1);
parms_constrained(2) = parms(2).^2;
parms_constrained(3:4) = parms(3:4).^2./(1+parms(3:4).^2);


