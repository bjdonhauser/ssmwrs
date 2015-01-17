function [h0] = calc_ss_h0(garchParms,p,q)
% calc_ss_h0    --  Calculates the steady-state variance of a GARCH(p,q)
%                   disturbance term.
%
%****f* SSMWRS/calc_ss_h0
%
% NAME
%   calc_ss_h0  --  Calculates the steady-state variance of a GARCH(p,q)
%                   disturbance term.
%
% SYNOPSIS
%   [h0] = calc_ss_h0(garchParms,p,q)
%
% INPUTS
%   * garchParms   --  [alpha_0, alpha_1, ..., alpha_p, alpha_(p+1), ...
%                       alpha_(p+q+1)] as described in the "description"
%                       below. (vector of length (p+q+1))
%                       Must only pass parmaters that produce STATIONARY
%                       variance.
%   * p             --  the "G" (Generalized) order of the GARCH
%                       disturbance. This corresponds to the h's below.
%   * q             --  the "AR" (Autoregressive) order of the GARCH
%                       distrubance. This corresponds to the AR
%
% OUTPUTS
%   * h0        --  steady-state variance of a GARCH(p,q)
%                   disturbance term.
%
% SIDE EFFECTS
%
% DESCRIPTION
%   If the disturbance term e_t ~ GARCH(p,q), then e_t ~ N(0,h_t) with
%   h_t = alpha_0 + alpha_1*h_(t-1) + ... + alpha_p*h_(t-p)
%                   alpha_(p+1)*e_(t-1)^2 + ... + alpha_(p+q+1)*e_(t-q)^2
%
%   NOTE: this is a little counterintuitive since h_t appears on the LHS of
%   the equation. So, we would think that this is the AR part when it is in
%   fact the G part.
%         
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   kalman_filter_garch
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-13
%
%***

nGarchParms = length(garchParms);
if nGarchParms ~= (p+q+1)
    error('GARCH(p,q) order does not match the number of passed parameters');
end

h0 = garchParms(1)/(1-sum(garchParms(2:nGarchParms)));



