function SSF_star = ucgarch11_ssf(parms, y)
% ucgarch11_ssf -- Takes the parameters [alpha_0, alpha_1, alpha_2,
%                  gamma_0, gamma_1, gamma_2] along with the data y and
%                  places it into the augmented state space form as seen in
%                  Kim and Nelson (1999, pg. 142).
%
%****f* SSMWRS/models/ucgarch11/ucgarch11_ssf
%
% NAME
%   ucgarch11_ssf -- Takes the parameters [alpha_0, alpha_1, alpha_2,
%                    gamma_0, gamma_1, gamma_2] along with the data y and
%                    places it into the augmented state space form as 
%                    seen in Kim and Nelson (1999, pg. 142).
%
% SYNOPSIS
%   [SSF_star] = ucgarch11_ssf(parms, y)
%
% INPUTS
%   * parms         -- corresponds to the values [alpha_0, alpha_1, 
%                      alpha_2, gamma_0, gamma_1, gamma_2]. (vector of
%                      length 6)
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
%   * SSF_star           -- a state space form object
%
% SIDE EFFECTS
%
% DESCRIPTION
%   This setup follows the augmented SSF_star of Kim and Nelson (1999, pg. 142)
%   and NOT the SSF with GARCH on pg. 140. 
%   FUTURE VERSION: We will just specify the SSF with GARCH and have a
%   program convert that to an augmented SSF object.
%
% REFERENCES
%   (1) Kim, C.-J. and Nelson, C. R., (1999), "State-Space Models with Regime
%   Switching", London: The MIT Press
%
% SEE ALSO
%   ucgarch11_main, ucgarch11_constraints, kalman_filter_constraints
%
% AUTHOR
%   Brian Donhauser
%
% CREATION DATE
%   2007-08-13
%
%***

k = 1; %length of the state vector
n = 1; %number of measurement equations

%=========================================================================%
% Unpack Hyperparameters
%=========================================================================%
alpha_0 = parms(1);
alpha_1 = parms(2);
alpha_2 = parms(3);
gamma_0 = parms(4);
gamma_1 = parms(5);
gamma_2 = parms(6);


%=========================================================================%
% Set up the state-space form on pg. 140
%=========================================================================%
SSF.H = 1;
SSF.A = 0;
SSF.R = 0;
SSF.mu = 0;
SSF.F = 1;
SSF.Q = 0;
SSF.LAMBDA = 1; % denotes the lambda in (6.4)
SSF.lambda = 1; % denotes the lambda in (6.5)

%=========================================================================%
% Starting Values
%=========================================================================%
SSF_star.h_me_0 = calc_ss_h0([alpha_0,alpha_1,alpha_2]',1,1);
SSF_star.h_te_0 = calc_ss_h0([gamma_0,gamma_1,gamma_2]',1,1);

SSF_star.b_0 = zeros((k+2),1);

% The state variable b_t = [tau_t, epsilon_t, omega_t]' where tau_t is the
% I(1) trend and epsilon_t and omega_t are GARCH disturbances with
% steady-state variance given by h10, h20 respectively. Thus, we enter the
% values for P0 as follows
SSF_star.P_0 = [eye(k)*1e+7, zeros(k,1), zeros(k,1); ...
               zeros(1,k), SSF_star.h_me_0, 0; ...
               zeros(1,k), 0, SSF_star.h_te_0];

%=========================================================================%
% Set up the state-space form on pg. 142
%=========================================================================%
SSF_star.F = [SSF.F, zeros(k,1), zeros(k,1); ...
              zeros(1,k),0,0; ...
              zeros(1,k),0,0];
SSF_star.mu = [SSF.mu; 0; 0];
SSF_star.G = [eye(k), zeros(k,1), SSF.lambda; ...
              zeros(1,k), 1, 0; ...
              zeros(1,k), 0, 1];
SSF_star.Q = [SSF.Q, zeros(k,1), zeros(k,1); ...
              zeros(1,k), SSF_star.h_me_0, 0; ...
              zeros(1,k), 0, SSF_star.h_te_0];  %NOTE: the h10 and h20 are
                                                %immaterial here since Q
                                                %will be updated at the
                                                %first iteration of the
                                                %Kalman Filter
SSF_star.H = [SSF.H, SSF.LAMBDA, zeros(n,1)];

%The others are the same
SSF_star.A = SSF.A;
SSF_star.R = SSF.R;


%=========================================================================%
% Data
%=========================================================================%
SSF_star.y = y;
SSF_star.z = zeros(length(y),1);

%=========================================================================%
% Parms
%=========================================================================%
% Must also pack these because they will need to be accessed in the h
% approximation step
SSF_star.parms = parms;






