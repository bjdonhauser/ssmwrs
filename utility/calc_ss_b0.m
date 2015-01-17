function [b0] = calc_ss_b0(F, mu)
% calc_ss_b -- Performs the calculation in equation (3.22) in Kim and
%              Nelson (1999, pg. 27). This function should only be called
%              if the state vector is indeed stationary and, hence, there
%              is actually a steady-state value for the state vector.
%
%****f* SSMWRS/calc_ss_b0
%
% NAME
%   calc_ss_b0  --  Performs the calculation in equation (3.22) in Kim and
%                   Nelson (1999, pg. 27). 
%
% SYNOPSIS
%   [b0] = calc_ss_b(F, mu)
%
% INPUTS
%   * F         --  The state transition matrix as indicated in equation 
%                   (3.2) in Kim and Nelson (1999, pg. 19). (Matrix of
%                   dimensions lStateVecXlStateVec
%   * mu        --  The state mean vector as indicated in equation (3.2) in
%                   Kim and Nelson (1999, pg. 19). (Vector of length 
%                   lStateVec)
%
% OUTPUTS
%   * b0        --  The steady-state value of the state-vector (vector of
%                   length lStateVec).
%
% SIDE EFFECTS
%
% DESCRIPTION
%   Performs the calculation in equation (3.22) in Kim and Nelson (1999, 
%   pg. 27). This function should only be called if the state vector is 
%   indeed stationary and, hence, there is actually a steady-state value 
%   for the state vector.
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

lStateVec = length(mu);

b0 = (eye(lStateVec)-F)\mu;
