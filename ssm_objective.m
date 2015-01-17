function [negLogL, results] = ssm_objective(parms, y, model_ssf, ...
                                               model_constraints, logL_start, filter)
% ssm_objective -- The objective function that the function
%                     ssm_optimizer optimizes. 
%
%****f* SSMWRS/ssm_objective
%
% NAME
%   ssm_objective -- The objective function that the function
%                       ssm_optimizer optimizes. 
%
% SYNOPSIS
%   [results] = ssm_objective(y, model_ssf, model_constraints, 
%                           kim_filter_0, parms_0, z, opt_engine)
%
% INPUTS
%   * parms             -- specifies the hyper paramters of the MS-SSM.
%                          (vector of length nParms)
%                          LATER RELEASE: will allow this to be an optional
%                          argument specifying initial paramters of zero.
%   * y                 -- the measurement variable. (vector length N)
%   * model_ssf         -- specifies the function which takes the
%                          hyperparamters of a model and returns a SSF 
%                          structure. (string) 
%   * model_constraints -- specifies the function which takes the
%                          hyperparamters and constraints them. If no
%                          constraints are necessary, a degenerate
%                          constraint function should be specified.
%                          (string)
%                          LATER RELEASE: will allow this to be an 
%                          optional argument.
%   * START             -- specifies value at which the logL function
%                          begins being recorded (scalar integer greater 
%                          than or equal to one)
%                          DEFAULT: 1 (i.e., the whole sample contributes
%                                   to the logL function)
%   * z                 -- FOR LATER RELEASE!!!
%                          auxilliary data. Refer to the state-space form
%                          used by Kim and Nelson (1999, pg. 100). (vector
%                          of length N).
%                          LATER RELEASE: will allow this to be an optional
%                          argument, with it defaulting to a vector of
%                          zeros.
%   * kim_filter_0      -- FOR LATE RELEASE!!!
%                          specifies the starting values for (1) the state
%                          variable, (2) the covariance of the state
%                          variable, and (3) the state probabilities.
%                          (structure).
%
% OUTPUTS
%   * negLogL           -- the negative of the logL value. The optimizer
%                          minimizes this value in order to maximize the
%                          log likelihood function.
%   * results           -- a strucure of results
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
%   
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   ssm_optimizer, kim_filter, lam_ssf, lam_constraints, lam_main
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-11
%
%***

% Constrain parameters

parms_constrained = model_constraints(parms);

% Aquire the SSF
% NOTE: will have to update model_ssf in order to take in kim_filter_0
SSF = model_ssf(parms_constrained, y);

% Call the kim filter
results = filter(SSF, logL_start);
negLogL = -results.logL;




