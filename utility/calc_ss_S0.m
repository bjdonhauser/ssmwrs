function [S0] = calc_ss_S0(p)
% calc_ss_S0    --  Performs the calculation in equation (4.49) in Kim and
%                   Nelson (1999, pg. 71). 
%
%****f* SSMWRS/calc_ss_S0
%
% NAME
%   calc_ss_S0  --  Performs the calculation in equation (3.26) in Kim and
%                   Nelson (1999, pg. 27). 
%
% SYNOPSIS
%   [S0] = calc_ss_S0(P)
%
% INPUTS
%   * p         --  The transition probability matrix as indicated in 
%                   equation (4.44) in Kim and Nelson (1999, pg. 71). 
%                   (Matrix of dimensions nStateXnState)
%
% OUTPUTS
%   * S0        --  The steady-state transition probabilities vector (a
%                   1XnState cell array of scalars)
%
% SIDE EFFECTS
%
% DESCRIPTION
%   Performs the calculation in equation (4.49) in Kim and Nelson (1999,
%   pg. 71).
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

nState = length(p);
A = [eye(nState) - p; ones(1, nState)];
S0 = (A'*A)\A';
S0 = S0(:,(nState+1))';
S0 = mat2cell(S0,1,ones(1,nState));



