function [results] = ssm_optimizer(parms_0, y, model_ssf, ...
                                    model_constraints, opt_engine, logL_start, opt, filter)
                      
% ssm_optimizer -- The main engine used to estimate markov-switching
%                     state-space models. This is intended to be run from a
%                     file 'model_main' where "model" is the identifier for
%                     a particular model (e.g., "lam"). 
%
%****f* SSMWRS/ssm_optimizer
%
% NAME
%   ssm_optimizer -- The main engine used to estimate markov-switching
%                       state-space models.
%
% SYNOPSIS
%   [results] = ssm_optimizer(y, model_ssf, model_constraints, 
%                           kim_filter_0, parms_0, z, opt_engine)
%
% INPUTS
%   * parms_0           -- specifies starting values for the paramters.
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
%                          LATER RELEASE: will allow this to be an optional
%                          argument, with it defaulting to a vector of
%                          zeros.
%   * opt_engine        -- specifies the optimization engine used (e.g., 
%                          fminsearch() and fminunc()). (string).
%                          LATER RELEASE: will allow this to be optional,
%                          with it defaulting to fminsearch()
%   * logL_start        -- specifies value at which the logL function
%                          begins being recorded (scalar integer greater 
%                          than or equal to one)
%                          DEFAULT: 1 (i.e., the whole sample contributes
%                                   to the logL function)
%   * opt               -- optimizer options.
%   * z                 -- FOR LATER RELEASE!!!
%                          auxilliary data. Refer to the state-space form
%                          used by Kim and Nelson (1999, pg. 100). (vector
%                          of length N). 
%   * kim_filter_0      -- FOR LATER RELEASE!!!
%                          specifies the starting values for (1) the state
%                          variable, (2) the covariance of the state
%                          variable, and (3) the state probabilities.
%                          (structure).
%
% OUTPUTS
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
%   ssm_objective, kim_filter, lam_ssf, lam_constraints, lam_main
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-11
%
%***

%========================================================================%
% Set defaults
%========================================================================%

if nargin <= 5 || isempty(logL_start)
    logL_start=1;
end

%========================================================================%
% Call the optimizer
%========================================================================%
tic;
[parms_c, fval, exitflag, output] = opt_engine('ssm_objective', ...
                                  parms_0, opt, y, model_ssf, ...
                                  model_constraints, logL_start, filter);
timer = toc;

%========================================================================%
% Calculate the standard errors and t-values
%========================================================================%
hess0 = fdhess('ssm_objective',parms_c, y, model_ssf, ...
                model_constraints, logL_start, filter);
cov0 = inv(hess0);
grad0 = gradnt(model_constraints, parms_c, 0.001);
cov1 = grad0*cov0*grad0';
parms_se = sqrt(diag(cov1));

%========================================================================%
% Recover the parameter values
%========================================================================%
parms_est = model_constraints(parms_c);

%========================================================================%
% Recover filter output and collect output
%========================================================================%
SSF = model_ssf(parms_est, y);
results = filter(SSF, logL_start)
results.parms_est = parms_est;
results.parms_se = parms_se;
results.exitflag = exitflag;
results.output = output;




