function [kim_filter_out] = kim_filter(SSF, logL_start)
% kim_filter -- 
%
%****f* SSMWRS/kim_filter
%
% NAME
%   kim_filter --  
%
% SYNOPSIS
%   [SSF] = kim_filter(SSF, [logL_start])
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
%   * kim_filter_out        -- a structure with the following fields:
%                  * logL   -- the log-likelihood value
%                  * S_t_j  -- Pr{S_t=j|I_t}.
%                              holds the filtered state inferences
%                              ((1XnObs) cell-array of 1xnStates cell
%                              arays.)
%                  * S_l_j  -- Pr{S_t=j|I_(t-1)}.
%                  * b_t    -- E{b_t|I_t}
%                              = Sum_j( E{b_t|S_t=j,I_t} *
%                                Pr{S_t=j|I_t} )
%                              = Sum_j( b_tt_j * S_t_j )
%                  * b_l    -- E{b_t|I_(t-1)}
%                              = Sum_j( E{b_t|S_t=j,I_(t-1)} *
%                                Pr{S_t=j|I_(t-1)} )
%                              = Sum_j( b_tt_j * S_l_j )
%                  * P_t    -- Cov{b_t|I_t}
%                  * P_l    -- Cov{b_t|I_(t-1)}
%                  * y_l    -- Forecast of y_t given info up to t-1
%                              = E{y_t|I_(t-1)}
%                              = Sum_j( ( H_j*E{b_t|I_(t-1),S_t=j} +
%                              A_j*z_t ) * Pr{S_t=j|I_(t-1)} )
%                              = DoubleSum_ij( y_tl_ij * 
%                                Pr{S_t=j,S_(t-1)=i|I_(t-1)} )
%                  * eta_l    -- prediction error
%                              = y_t - y_tl
%                  * f_l    -- conditional variance of prediction error
%                              = E{eta_tl*eta_tl'|I_(t-1)}
%                              = DoubleSum_ij( f_tl_ij * 
%                                Pr{S_t=j,S_(t-1)=i|I_(t-1)} )
%   * logL_start            -- (optional) indicates how many data points in
%                              before the log-likelihood value begins being
%                              recorded. This may be added if we don't want
%                              the initial forecast errors contributing the
%                              log-likelihood function.
%                              Default: 1 (indicating the entire sample)
%                              LATER RELEASE: determine whether the Kim
%                              Filter should allow specification of this or
%                              not.
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
%   p = [p11, p21, ... , pM1; p12, ..., pM2; ... ; p1M, ... pMM]; where pij
%   = Pr{S_t=j|S_(t-1)=i}
%
%   The algrorithm follows the diagram of Kim and Nelson (1999, pg. 105)
%   and the description in Kim and Nelson (1999, pg. 104). 
%
%
%
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-08
%
%***

%========================================================================%
% Important Variables
%========================================================================%

%   * b_tt_j -- E{b_t|S_t=j,I_t}                     (5.8) 
%               holds the filtered state variable inferences
%               ((1XnObs) cell-array of (1XnStates) cell
%               arrays of (lStateVectorX1) vectors).
%   * b_tl_i -- E{b_t|I_(t-1),S_(t-1)=i}
%               = Sum_j( E{b_t|I_(t-1),S_(t-1)=i,S_t=j} *
%                 Pr{S_t=j|I_(t-1)} )
%               = Sum_j( b_tl_ij * S_l_j )
%   * P_tt_j -- Cov{b_t|S_t=j,I_(t-1)}               (5.9) 
%               holds the filtered state variable covariance
%               matrices ((1XnObs) cell-array of (1XnStates) cell
%               arrays of (lStateVectorXlStateVector)
%               matrices).
%   * P_tl_i -- Cov{b_t|I_(t-1),S_(t-1)=i}

%========================================================================%
% Set Defaults
%========================================================================%

if nargin <= 1
    logL_start = 1;
end

%========================================================================%
% Variable lengths
%========================================================================%
nObs = length(SSF.y);
lStateVec = length(SSF.b0{1});
nStates = length(SSF.p);
nMsrmntEqns = length(SSF.R{1});

%========================================================================%
% Initialize
%========================================================================%
% NOTE: nObs+1 to take into account the intitial starting value
logL = 0;
S_t_j = cell(1,nObs+1);
S_l_j = cell(1,nObs);

b_l = cell(1,nObs);
P_l = cell(1,nObs);
y_l = cell(1,nObs);
eta_l = cell(1,nObs);
f_l = cell(1,nObs);

b_t = cell(1,nObs);
P_t = cell(1,nObs);

%========================================================================%
% Starting Values
%========================================================================%
S_t_j{1} = SSF.S0;
b_tt_j = SSF.b0;        
P_tt_j = SSF.P0;         


for t=1:nObs
    
    %=====================================================================%
    % Kalman Filter (equations 5.10-5.15, pg. 100 (plus an extra step) )
    %=====================================================================%
    for i=1:nStates
        for j=1:nStates
           
            b_tl_ij{i}{j} = SSF.mu{j} + SSF.F{j}*b_tt_j{i};
            
            P_tl_ij{i}{j} = SSF.F{j}*P_tt_j{i}*SSF.F{j}' +  ...
                            SSF.G{j}*SSF.Q{j}*SSF.G{j}';
            
            % EXTRA: to get the forecast            
            y_tl_ij{i}{j} = SSF.H{j}*b_tl_ij{i}{j} + SSF.A{j}*SSF.z(t); 
            
            % NOTE: different from book because of the extra step
            eta_tl_ij{i}{j} = SSF.y(t) - y_tl_ij{i}{j};
            
            f_tl_ij{i}{j} = SSF.H{j}*P_tl_ij{i}{j}*SSF.H{j}'+SSF.R{j};
            
            b_tt_ij{i}{j} = b_tl_ij{i}{j} + P_tl_ij{i}{j}*SSF.H{j}'* ...
                            (f_tl_ij{i}{j}\eta_tl_ij{i}{j});
                        
            P_tt_ij{i}{j} = (eye(lStateVec) - P_tl_ij{i}{j}*SSF.H{j}'* ...
                            (f_tl_ij{i}{j}\SSF.H{j}))*P_tl_ij{i}{j};
            
        end
    end

    %==================================================================%
    % Hamilton Filter (equations 5.18-5.23 in Kim and Nelson (1999, pgs.
    % 102-103) 
    %==================================================================%
    fy_l = 0;
    for i=1:nStates
        for j=1:nStates
            
            %============================================================%
            % Steps 1 and 2 (Kim and Nelson (1999, pg. 102))
            %============================================================%       
            % NOTE: p(j,i) and not p(i,j). See equation (5.4) for notation details.
            % (5.18) Pr{S_t=j,S_(t-1)=i|I_(t-1)} = Pr{S_t=j|S_(t-1)=i}*Pr{S_(t-1)=i|I_(t-1)}
            S_l_ij{i}{j} = SSF.p(j,i)*S_t_j{t}{i};   
            
            %(5.21) f(y_t|S_t=j,S_(t-1)=i,I_(t-1)) = ...
            fy_l_ij_c = ((2*pi)^(-nMsrmntEqns/2)) * ((det(f_tl_ij{i}{j}))^(-0.5)) * ...
                        exp(-0.5*eta_tl_ij{i}{j}'* ...
                        (f_tl_ij{i}{j}\eta_tl_ij{i}{j}));
            
            % f(y_t,S_t=j,S_(t-1)=i|I_(t-1)) =
            % f(y_t|S_t=j,S_(t-1)=i,I_(t-1)) * Pr{S_t=j,S_(t-1)=i|I_(t-1)}
            fy_l_ij_uc{i}{j} = fy_l_ij_c * S_l_ij{i}{j};
            
            % (5.20) f(y_t|I_(t-1)) = DoubleSum( f(y_t,S_t=j,S_(t-1)=i|I_(t-1)) )           
            fy_l = fy_l + fy_l_ij_uc{i}{j};
            %============================================================%
            
        end       
    end
    
    % (5.24) log(theta) = log(theta) + log(f(y_t|I_(t-1))
    if t >= logL_start
        logL = logL + log(fy_l);
    end
    
    %==============================================%
    % Step 3 (Kim and Nelson (1999, pg. 103))
    %==============================================%
    % NOTE: There are extra steps taken here to get Pr{S_t=j|I_(t-1)} in
    % addition to Pr{S_t=j|I_t}
    
    for i=1:nStates
        for j=1:nStates            
            % (5.22) Pr{S_t=j,S_(t-1)=i|I_t} =
            % f(y_t,S_t=j,S_(t-1)=i|I_(t-1)) / f(y_t|I_(t-1))
            S_t_ij{i}{j} = fy_l_ij_uc{i}{j} / fy_l;         
        end
    end
    
    % NOTE: the reversal of the order of iteration, which is necessary.
    for j=1:nStates        
        
        S_t_j{t+1}{j} = 0;      % NOTE: first entry not empty because of 
                                % the initial values
        S_l_j{t+1}{j} = 0;      % NOTE: the first entry will be empty
        
        for i=1:nStates            
            
            % (5.23) Pr{S_t=j|I_t} = sum_i( Pr{S_(t-1)=i,S_t=j|I_t} )          
            S_t_j{t+1}{j} = S_t_j{t+1}{j} + S_t_ij{i}{j};
            
            % Pr{S_t=j|I_(t-1)} = sum_i( Pr{S_(t-1)=i,S_t=j|I_t} ) 
            S_l_j{t+1}{j} = S_l_j{t+1}{j} + S_l_ij{i}{j};
        
        end        
    end
    
    %==================================================================%
    % Collapsing Terms (equations 5.16)
    %==================================================================%
    for j=1:nStates     
        b_tt_j{j} = zeros(lStateVec,1);
        for i=1:nStates
            % (5.16) b_j_tt = Sum( Pr{S_(t-1)=i,S_t=j|I_t}*b_tt_ij ) /
            % Pr{S_t=j|I_(t-1)}
            b_tt_j{j} = b_tt_j{j} + S_t_ij{i}{j} * b_tt_ij{i}{j};
        end
        b_tt_j{j} = b_tt_j{j} / S_t_j{t+1}{j};
    end

    %==================================================================%
    % Collapsing Terms (equations 5.17)
    %==================================================================%    
    for j=1:nStates
        P_tt_j{j} = zeros(lStateVec,lStateVec);
        for i=1:nStates
            % (5.17') P_j_tt = Sum( Pr{S_(t-1)=i,S_t=j|I_t} * {P_tt_ij +
            % (b_j_tt - b_tt_ij)*(b_j_tt - b_tt_ij)'} ) / Pr{S_t=j|I_t}
            P_tt_j{j} = P_tt_j{j} + S_t_ij{i}{j} * ( P_tt_ij{i}{j} + ...
                ( b_tt_j{j} - b_tt_ij{i}{j} ) * ( b_tt_j{j} - b_tt_ij{i}{j} )' );
         end
         P_tt_j{j} = P_tt_j{j}/S_t_j{t+1}{j};
    end
    
    %==================================================================%
    % Collecting output
    %==================================================================%
    
    b_l{t} = zeros(lStateVec,1);
    P_l{t} = zeros(lStateVec,lStateVec);
    y_l{t} = zeros(nMsrmntEqns,1);
    eta_l{t} = zeros(nMsrmntEqns,1);
    f_l{t} = zeros(nMsrmntEqns,nMsrmntEqns);
        
    b_t{t} = zeros(lStateVec,1);
    P_t{t} = zeros(lStateVec,lStateVec);
    
    for i=1:nStates       
        for j=1:nStates
            b_l{t} = b_l{t} + b_tl_ij{i}{j}*S_l_ij{i}{j};
            P_l{t} = P_l{t} + P_tl_ij{i}{j}*S_l_ij{i}{j};
            y_l{t} = y_l{t} + y_tl_ij{i}{j}*S_l_ij{i}{j};
            eta_l{t} = eta_l{t} + eta_tl_ij{i}{j}*S_l_ij{i}{j};
            f_l{t} = f_l{t} + f_tl_ij{i}{j}*S_l_ij{i}{j};

            b_t{t} = b_t{t} + b_tt_ij{i}{j}*S_t_ij{i}{j};
            P_t{t} = P_t{t} + P_tt_ij{i}{j}*S_t_ij{i}{j};           
        end        
    end
    
    
end

%============================================%
% Fill output structure
%============================================%
kim_filter_out.logL = logL;

kim_filter_out.S_t_j = S_t_j(2:(nObs+1));
kim_filter_out.S_l_j = S_l_j(2:(nObs+1));       %This is OK because of how I initialized

kim_filter_out.b_t = b_t;
kim_filter_out.b_l = b_l;

kim_filter_out.P_t = P_t;
kim_filter_out.P_l = P_l;

kim_filter_out.y_l = y_l;
kim_filter_out.eta_l = eta_l;
kim_filter_out.f_l = f_l;







