function SSF = lam_ssf(parms, y)
% lam_ssf -- Takes the 7 parameters from the Lam model of output and
% places them into state-space form as indicated in Kim and Nelson (1999,
% pgs. 111-112).
%
%****f* SSMWRS/lam_ssf
%
% NAME
%   lam_ssf -- Takes the 7 parameters from the Lam model of output and
%   places them into state-space form as indicated in Kim and Nelson (1999,
%   pgs. 111-112).
%
% SYNOPSIS
%   [SSF] = lam_ssf(parms, y)
%
% INPUTS
%   * parms         -- vector of length 7 corresponding to 
%                      [p, q, delta_0, delta_1, sigma, phi_1, phi_2] of 
%                      the Lam model 
%   * y             -- the data. Should be the pre-processed growth rate 
%                      of gdp (first difference of log(gdp) a vector of 
%                      length N.
%   * z             -- FOR LATER RELEASE!!!
%                      the auxiliary data. In this model, it's just a
%                      vector ones. (a vector of length N)
%   * kim_filter_0  -- FOR LATER RELEASE!!!
%                      specifies the starting values for (1) the state
%                      variable, (2) the covariance of the state
%                      variable, and (3) the state probabilities.
%                      (structure).
%
% OUTPUTS
%   * SSF           -- a state space form object
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
%   The Lam Model of output growth is represented (in abbreviated form) as:
%   the following in Kim and Nelson (1999, pgs. 111-112).
%
%   Measurement equation:
%   d_y_t = [1 1]*[x_t x_(t-1)]' + delta_j
%   where d_y_t, is the growth rate of gdp, x_t represents the stationary
%   cycle component of gdp, and delta represents the markov-switching mean
%   growth rate of gdp
%
%   Transition equation:
%   [x_t x_(t-1)]' = [theta_1 theta_2; 1, 0]*[x_(t-1) x_(t-2)]'+[u_t, 0]'
%   where [u_t 0]' ~ N([0, 0]',[s2, 0; 0, 0])
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   Clark_SSF.m (Author: Drew Creal),
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-10
%
%***

%=========================================================================%
% Unpack Hyperparameters
%=========================================================================%
p = parms(1);
q = parms(2);
delta_1 = parms(3);
delta_2 = parms(4);
sigma = parms(5);
phi_1 = parms(6);
phi_2 = parms(7);


%=========================================================================%
% Set up the state-space form with hyperparamters
%=========================================================================%
SSF.H{1} = [1, -1];
SSF.H{2} = [1, -1];

SSF.A{1} = [delta_1];  
SSF.A{2} = [delta_2];   

SSF.R{1} = [0];
SSF.R{2} = [0];

SSF.mu{1} = [0];
SSF.mu{2} = [0];

SSF.F{1} = [phi_1, phi_2; 1, 0];
SSF.F{2} = [phi_1, phi_2; 1, 0];

SSF.G{1} = eye(2);
SSF.G{2} = eye(2);

SSF.Q{1} = [sigma^2, 0; 0, 0];
SSF.Q{2} = [sigma^2, 0; 0, 0];

SSF.p = [q, (1-p); (1-q), p];

if (any(abs(sum(SSF.p)-1) > 1e-6))
    error('Columns must add to 1. See notation of Kim and Nelson (1999, pg. 98)');
end

%=========================================================================%
% Starting Values
%=========================================================================%

SSF.P0{1} = calc_ss_P0(SSF.F{1},SSF.Q{1});
SSF.P0{2} = SSF.P0{1};
SSF.S0 = calc_ss_S0(SSF.p);
SSF.b0{1} = zeros(2,1);
SSF.b0{2} = SSF.b0{1};

% % The Kim and Nelson text and Kim (1994) starting values. Different from
% % those used in kim_je.opt.
% SSF.b0{1} = [5.224, 0.535]';  
% SSF.b0{2} = [5.224, 0.535]'; 
% 
% % The kim_je.opt starting values;
% % SSF.b0{1} = zeros(2,1); 
% % SSF.b0{2} = SSF.b0{1}; 
% 
% % The Kim and Nelson (1999) and Kim (1994) starting values. Different 
% % from those used in kim_je.opt.
% SSF.P0{1} = zeros(2,2);
% SSF.P0{2} = SSF.P0{1};
% 
% % The kim_je.opt starting values.
% % SSF.P0{1} = [5.5541, 5.1741; 5.1741, 5.5541];
% % SSF.P0{2} = SSF.P0{1};
% 
% % Diffuse starting values
% %SSF.P0{1} = eye(2)*1e+7;  
% %SSF.P0{2} = eye(2)*1e+7;
% 
% SSF.S0{1} = (1-p)/(2-p-q); 
% SSF.S0{2} = 1-SSF.S0{1};

%=========================================================================%
% Data
%=========================================================================%
SSF.y = y;
SSF.z = ones(length(y),1);
