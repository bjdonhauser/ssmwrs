function [kalman_filter_out] = kalman_filter(SSF, logL_start)
% kalman_filter -- 
%
%****f* SSMWRS/kalman_filter
%
% NAME
%   kalman_filter --  
%
% SYNOPSIS
%   [SSF] = kalman_filter(SSF, [logL_start])
%
% INPUTS
%   * SSF       -- the state-space form structure (see _parms2ssf for
%                  details). This includes the data y and z and intial 
%                  values.
%   * logL_start-- (optional). Specifies iteration at which the logL begins
%                  being recorded.
%                  Default: 1. (i.e., the loglikelihood function is
%                  recorded for the entire sample)
%
% OUTPUTS
%   * kalman_filter_out     -- a structure with the following fields:
%       * logL                  -- the log-likelihood value
%       * b_t                   -- E{b_t|I_t}
%       * P_t                   -- Cov{b_t|I_t}
%       * b_l                   -- E{b_t|I_(t-1)}
%       * P_l                   -- Cov{b_t|I_(t-1)}
%       * y_l                   -- Forecast of y_t given info up to t-1
%       * h_l                   -- prediction error
%                                  = y_t - y_tl
%       * f_l                   -- conditional variance of prediction error
%                                  = E{h_tl*h_tl'|I_(t-1)}
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
%   A particular general form of a state-space model is given in Kim and
%   Nelson (1999, pg. 20-21, eqns (3.1)-3.4)). To stay consistent with 
%   chapter 5 of the same text, we'll instead adopt the general form of 
%   the model there (pg. 98, eqns (5.1)-(5.3)) and drop the state notation
%   "S_t". Thus, we get the following general form:
%
%   Measurement equation:
%   y_t = H*b_t + A*z_t + e_t           e_t ~ N(0, R)
%
%   Transition equation:
%   b_t = mu + F*b_(t-1) + G*v_t    v_t ~ N(0, Q)
%
%   The prediction and updating equations follow the notation of pg. 100,
%   eqns (5.10)-(5.15) while dropping the "(i,j)" notation. 
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
k       = length(SSF.b0);   % length of the state vector
n       = length(SSF.R);    % the number of measurement equations

%========================================================================%
% Initialize
%========================================================================%
logL = 0;
b_t = cell(1,nObs);
P_t = cell(1,nObs);
b_l = cell(1,nObs);
P_l = cell(1,nObs);
y_l = cell(1,nObs);
h_l = cell(1,nObs);
f_l = cell(1,nObs);

%========================================================================%
% Starting Values
%========================================================================%
b_tt = SSF.b0;        
P_tt = SSF.P0;         

%========================================================================%
% Iterate
%========================================================================%
for t=1:nObs  
    %=====================================================================%
    % Kalman Filter (equations 5.10-5.15, pg. 100 without i,j notation )
    %=====================================================================%
    b_tl = SSF.mu + SSF.F*b_tt;
            
    P_tl = SSF.F*P_tt*SSF.F' + SSF.G*SSF.Q*SSF.G';
            
    % EXTRA: to get the forecast            
    y_tl = SSF.H*b_tl + SSF.A*SSF.z(t); 
            
    % NOTE: different from book because of the extra step
    h_tl = SSF.y(t) - y_tl;
            
    f_tl = SSF.H*P_tl*SSF.H'+SSF.R;
            
    b_tt = b_tl + P_tl*SSF.H'*(f_tl\h_tl);
                        
    P_tt = (eye(k) - P_tl*SSF.H'*(f_tl\SSF.H))*P_tl;
    
    %=====================================================================%
    % Add to log likelihood function
    %=====================================================================%
    
    fy      = ((2*pi)^(-n/2))*((det(f_tl)^(-0.5))*
              exp(-0.5*eta_tl'*(f_tl\eta_tl));
              
    logL    = logL + log(fy);
    
    %==================================================================%
    % Collecting output
    %==================================================================%
    b_t{t} = b_tt;
    P_t{t} = P_tt;
    b_l{t} = b_tl;
    P_l{t} = P_tl;
    y_l{t} = y_tl;
    h_l{t} = h_tl;
    f_l{t} = f_tl;
end

%============================================%
% Fill output structure
%============================================%
kim_filter_out.logL = logL;
kim_filter_out.b_t = b_t;
kim_filter_out.P_t = P_t;
kim_filter_out.b_l = b_l;
kim_filter_out.P_l = P_l;
kim_filter_out.y_l = y_l;
kim_filter_out.h_l = h_l;
kim_filter_out.f_l = f_l;









