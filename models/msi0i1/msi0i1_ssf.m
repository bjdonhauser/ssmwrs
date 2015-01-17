function SSF = msi0i1_ssf(parms, y)
% msi0i1_ssf -- Takes the parameters sig2eta, sig2eps_j, p11, and p22 and 
% returns the state space form object SSF for the msi0i1 model. 
%
%****f* SSMWRS/models/msi0i1/msi0i1_ssf
%
% NAME
%   msi0i1_ssf -- Takes the parameters sig2eta, sig2eps_j, p11, and p22 and 
%   returns the state space form object SSF for the msi0i1 model. 
%
% SYNOPSIS
%   [SSF] = msi0i1_ssf(parms)
%
% INPUTS
%   * parms         -- vector of length 4 corresponding to the values
%                      [sig2eta, sig2eps_j, p11, and p22]
%   * y             -- the measurement data. (a vector length N)
%   * z             -- FOR LATER VERSION!!!
%                      the auxiliary data. In this model, it's just a
%                      vector zeros. (a vector of length N)
%   * kim_filter_0  -- FOR LATER VERSION!!!
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
%   The msi0i1 model is intended to be an alternative to the Stock and
%   Watson (2006) UC-SV model of inflation. In the msi0i1 model there is
%   Markov switching between inflation having an I(0) and I(1) component.
%   The model looks as follows:
%
%   pi(t) = (1)*tau(t) + eta            eta ~ N(0,sig2eta)
%   tau(t) = (1)*tau(t-1) + epsilon     epsilon ~ N(0,sig2eps_j)
%
%   sigma_epsilon_j^2 = sigma_epsilon_1^2 when S_t = 2 and 0 when S_t = 1
%
%   Pr{S_t=2|S_(t-1)=2}=p22, Pr{S_t=1|S_(t-1)=1}=p11
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
%   2007-08-07
%
%***

%=========================================================================%
% Unpack Hyperparameters
%=========================================================================%
sig2eta = parms(1);
sig2eps_2 = parms(2);
p11 = parms(3);
p22 = parms(4);

%=========================================================================%
% Set up the state-space form with hyperparamters
%=========================================================================%
SSF.H{1} = [1];
SSF.H{2} = [1];

SSF.A{1} = [0];
SSF.A{2} = [0];

SSF.R{1} = [sig2eta];
SSF.R{2} = [sig2eta];

SSF.mu{1} = [0];
SSF.mu{2} = [0];

SSF.F{1} = [1];
SSF.F{2} = [1];

SSF.G{1} = [1];
SSF.G{2} = [1];

SSF.Q{1} = [0];
SSF.Q{2} = [sig2eps_2];

% NOTE: lowercase "p" corresponds to the transition matrix while upper case
% "P" will correspond to the mean-squared error forecast.
SSF.p = [p11, (1-p22); (1-p11), p22];

if (any(abs(sum(SSF.p)-1) > 1e-6))
    error('Columns must add to 1. See notation of Kim and Nelson (1999, pg. 98)');
end

%=========================================================================%
% Starting Values
%=========================================================================%
SSF.S0 = calc_ss_S0(SSF.p);

% % state1 is the stationary state and state2 is the non-stationary state.
% % The setting of the initial values roughly follows Kim and Nelson (1999,
% % pg. 27). In the stationary case we set the value of the state vector to
% % the intial value of the series. This is because the trend component acts
% % like a mean in this case. Also, the initial value of the state vector is
% % set to the initial observation in the non-stationary case
SSF.b0{1} = [0];   
SSF.b0{2} = [0];

% !!!NOTE: The state variable is just a mean in the stationary, state1
% case. This means that the covariance matrix, P, will just be the
% degenerate 0. This could be BIG TROUBLE. We'll set it here to be very
% small.

% !!!NOTE: This won't run. Degenerate case leads to division by zero!
%SSF.P0{1} = calc_ss_P0(SSF.F{1},SSF.Q{1});

SSF.P0{1} = 1e+7;
SSF.P0{2} = 1e+7;

%=========================================================================%
% Data
%=========================================================================%
SSF.y = y;
SSF.z = zeros(length(y),1);






