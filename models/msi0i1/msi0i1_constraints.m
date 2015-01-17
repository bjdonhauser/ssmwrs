function parms_constrained = msi0i1_constraints(parms)
% msi0i1_constraints -- Transforms the unconstrained
% parameters of the msi0i1 model to constrained parameters. The purpose of
% doing this is so that and unconstrained optimization can be performed in
% the new, transformed variables
%
%****f* SSMWRS/models/msi0i1/msi0i1_constraints
%
% NAME
%   msi0i1_constraints --  Transforms the unconstrained
% parameters of the msi0i1 model to constrained parameters. The purpose of
% doing this is so that and unconstrained optimization can be performed in
% the new, transformed variables  
%
% SYNOPSIS
%   [SSF] = msi0i1_constraints(parms)
%
% INPUTS
%   * parms             -- vector of length 4 corresponding to the values
%                          [sig2eta, sig2eps_j, p11, and p22]
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
%   The constraints are to keep sig2eps and sig2eta positive and p and q
%   positive and less than 1. 
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   ssm_optimizer, ssm_objective, msi0i1_main, msi0i1_ssf
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
parms_constrained(1:2) = parms(1:2).^2;
parms_constrained(3:4) = parms(3:4).^2./(1+parms(3:4).^2);


