function [kalman_filter_out] = kalman_filter_garch(SSF, logL_start)
% kalman_filter_garch -- Runs the kalman_filter on the augmented GARCH SSF
%                        object while doing the h approximations explained
%                        in Kim and Nelson (1999, pg. 143). 
%                        CURRENT: ONLY PROCESS GARCH(1,1)!!! for ME and TE
%
%****f* SSMWRS/kalman_filter_garch
%
% NAME
%   kalman_filter_garch --  Runs the kalman_filter on the augmented GARCH 
%                           SSF object while doing the h approximations 
%                           explained in Kim and Nelson (1999, pg. 143).
%
% SYNOPSIS
%   [SSF] = kalman_filter_garch(SSF, [logL_start])
%
% INPUTS
%   * SSF       -- the state-space form structure. This is a SPECIAL 
%                  state space form structure in that it follows Kim and 
%                  Nelson (1999, pg. 142, eqns (6.24'), (6.25), (6.26'), 
%                  and (6.27). This includes the data y and z and intial 
%                  values.
%   * logL_start-- (optional). Specifies iteration at which the logL begins
%                  being recorded.
%                  Default: 1. (i.e., the loglikelihood function is
%                  recorded for the entire sample)
%
% OUTPUTS
%   * kalman_filter_garch_out   -- a structure with the following fields:
%       * logL                  -- the log-likelihood value
%       * b_t                   -- E{b_t|I_t}
%       * P_t                   -- Cov{b_t|I_t}
%       * b_l                   -- E{b_t|I_(t-1)}
%       * P_l                   -- Cov{b_t|I_(t-1)}
%       * y_l                   -- Forecast of y_t given info up to t-1
%       * eta_l                   -- prediction error
%                                  = y_t - y_tl
%       * f_l                   -- conditional variance of prediction error
%                                  = E{eta_tl*eta_tl'|I_(t-1)}
%   * logL_start            -- (optional) indicates how many data points in
%                              before the log-likelihood value begins being
%                              recorded. This may be added if we don't want
%                              the initial forecast errors contributing the
%                              log-likelihood function.
%                              Default: 1 (indicating the entire sample)
%                              LATER RELEASE: determine whether the kalman
%                              Filter should allow specification of this or
%                              not.
%
% SIDE EFFECTS
%
% DESCRIPTION
%   Takes an augmented GARCH SSF object outlined in equations (6.24'), 
%   (6.25), (6.26'), and (6.27) in Kim and Nelson (1999, pg. 143) and runs 
%   the Kalman recursions of (6.28)-(6.33) while doing the "h"
%   approximations described thereafter.
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   kim_filter
%   
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-13
%
%***

%========================================================================%
% Set Defaults
%========================================================================%

if nargin <= 1
    logL_start = 1;
end

%========================================================================%
% Variable lengths
%========================================================================%
nObs    = length(SSF.y);
k       = length(SSF.b_0);      % NOTE: this is the length of the 
                                % augmented state vector. The original 
                                % state vector is of length k-2
n       = length(SSF.R);        % number of measurement equations

%========================================================================%
% Initialize
%========================================================================%
logL    = 0;
b_t     = cell(1,nObs);
P_t     = cell(1,nObs);
b_l     = cell(1,nObs);
P_l     = cell(1,nObs);
y_l     = cell(1,nObs);
eta_l   = cell(1,nObs);
f_l     = cell(1,nObs);
h_me_t  = cell(1,nObs);
h_te_t  = cell(1,nObs);

%========================================================================%
% Starting Values
%========================================================================%
b_tt = SSF.b_0;        
P_tt = SSF.P_0; 
h_me = SSF.h_me_0; 
h_te = SSF.h_te_0; 

%========================================================================%
% Iterate
%========================================================================%
for t=1:nObs  
    %=====================================================================%
    % Kalman Filter (equations 6.28-6.33, pg. 143)
    %=====================================================================%
    b_tl = SSF.mu + SSF.F*b_tt;
    
    %=====================================================================%
    % h appx. step explained on pg. 143
    %=====================================================================%            
    epsilon2_l  = b_tt(k-1)^2 + P_tt((k-1),(k-1));
    omega2_l    = b_tt(k)^2 + P_tt(k,k);
    
    h_me = SSF.parms(1) + SSF.parms(2)*epsilon2_l + ...
           SSF.parms(3)*h_me;
    h_te = SSF.parms(4) + SSF.parms(5)*omega2_l + ...
           SSF.parms(6)*h_te; 
    
    SSF.Q((k-1),(k-1))  = h_me;
    SSF.Q(k,k)          = h_te;
    
    %=====================================================================%
    % Resume Kalman Filter (equations 6.28-6.33, pg. 143)
    %=====================================================================%    
    P_tl    = SSF.F*P_tt*SSF.F' + SSF.G*SSF.Q*SSF.G';
            
    % EXTRA: to get the forecast            
    y_tl    = SSF.H*b_tl + SSF.A*SSF.z(t); 
            
    % NOTE: different from book because of the extra step
    eta_tl  = SSF.y(t) - y_tl;
            
    f_tl    = SSF.H*P_tl*SSF.H'+SSF.R;
            
    b_tt    = b_tl + P_tl*SSF.H'*(f_tl\eta_tl);
                        
    P_tt    = (eye(k) - P_tl*SSF.H'*(f_tl\SSF.H))*P_tl;
    
    %=====================================================================%
    % Add to log likelihood function
    %=====================================================================%
    
    fy      = ((2*pi)^(-n/2))*((det(f_tl)^(-0.5))* ...
              exp(-0.5*eta_tl'*(f_tl\eta_tl)));
              
    logL    = logL + log(fy);
    
    %==================================================================%
    % Collecting output
    %==================================================================%
    b_t{t}      = b_tt;
    P_t{t}      = P_tt;
    b_l{t}      = b_tl;
    P_l{t}      = P_tl;
    y_l{t}      = y_tl;
    eta_l{t}    = eta_tl;
    f_l{t}      = f_tl;
    h_me_t{t}   = h_me;
    h_te_t{t}   = h_te;
end

%============================================%
% Fill output structure
%============================================%
kalman_filter_out.logL     = logL;
kalman_filter_out.b_t      = b_t;
kalman_filter_out.P_t      = P_t;
kalman_filter_out.b_l      = b_l;
kalman_filter_out.P_l      = P_l;
kalman_filter_out.y_l      = y_l;
kalman_filter_out.eta_l    = eta_l;
kalman_filter_out.f_l      = f_l;
kalman_filter_out.h_me_t   = h_me_t;   %this is the ME GARCH distrubance variance
kalman_filter_out.h_te_t   = h_te_t;   %this is the TE GARCH distrubance variance