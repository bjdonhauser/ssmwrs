function [P0] = calc_ss_P0(F, Q)
% calc_ss_P0    --  Performs the calculation in equation (3.26) in Kim and
%                   Nelson (1999, pg. 27). This function should only be 
%                   called if the state vector is indeed stationary and, 
%                   hence, there is actually a steady-state value for the 
%                   covariance matrix.
%
%****f* SSMWRS/calc_ss_P0
%
% NAME
%   calc_ss_P0  --  Performs the calculation in equation (3.26) in Kim and
%                   Nelson (1999, pg. 27). 
%
% SYNOPSIS
%   [P0] = calc_ss_P0(F, mu)
%
% INPUTS
%   * F         --  The state transition matrix as indicated in equation 
%                   (3.2) in Kim and Nelson (1999, pg. 19). (Matrix of
%                   dimensions lStateVecXlStateVec)
%   * Q         --  The state error matrix as indicated in equation (3.4) 
%                   in Kim and Nelson (1999, pg. 19). (Matrix of
%                   dimensions lStateVecXlStateVec)
%
% OUTPUTS
%   * P0        --  The steady-state value of the covariance matrix of the
%                   state vector (vector of length lStateVec).
%
% SIDE EFFECTS
%
% DESCRIPTION
%   Performs the calculation in equation (3.26) in Kim and Nelson (1999, 
%   pg. 27). This function should only be called if the state vector is 
%   indeed stationary and, hence, there is actually a steady-state value 
%   for the covariance matrix.
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
%   2007-08-11
%
%***

lStateVec = length(F);

P0 = (eye(lStateVec^2)-kron(F,F))\reshape(Q,lStateVec^2,1);

P0 = reshape(P0,lStateVec,lStateVec);


