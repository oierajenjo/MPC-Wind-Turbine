clc
clear all
close all

%% Obtain all variables
variables_IPC
load('BladedFiles\performancemap_data.mat')
Constant_variables
addpath('functions');

% kAns = questdlg('Execute Kalman Filter with Bladed measurements', ...
% 	'Kalman Filter');
kAns = 'Yes';
kAns = convertCharsToStrings(kAns);
if kAns=="Cancel"
    return
end

u_b = ones(4,N);
theta_f = 0;
[lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
u_b = [theta_f; theta_f; theta_f; K].*u_b;


[f,h,Q,R] = equationIPC(var,Ts,ct_l,cp_l,lambdaVec,pitchVec,Lk);

% Step 3: Initialize state and covariance
% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(x) sqrt(Q(x))*randn(Lk, 1); % Generate random process noise (from assumed Q)
% v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
v = zeros(Yk, N);

%% Runge-Kutta 4th order method
% Initialize matrices
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k
[xt,yt] = RK4(f,xt,u_b,h,yt,n,v,var,Ts);

true_plots(yt,y_me,xt,data,t)

%% Unscented Kalman Filter
% Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i;
P0 = [M.sigma_enc; M.sigma_tdef; M.sigma_tvel; M.sigma_tdef; M.sigma_tvel;...
    M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel; M.sigma_bvel;...
    M.sigma_bvel; M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel;...
    M.sigma_bvel; M.sigma_bvel; M.sigma_pit; M.sigma_pit; M.sigma_pit;...
    M.sigma_pitvel; M.sigma_pitvel; M.sigma_pitvel; M.sigma_pow;...
    M.sigma_vane; 0.01; M.sigma_azim].^2;
P0 = diag(P0);
% P0 = 0.01*eye(Lk,Lk);

if kAns == "Yes"
    disp("Measured values")
    [xk,P,e] = UKF(f,h,Q,R,xk,y_me,u_b,Lk,Yk,P0,Ts,v,n);
elseif kAns == "No"
    disp("True values")
    [xk,P,e] = UKF(f,h,Q,R,xk,yt,u_b,Lk,Yk,P0,Ts,v,n);
end
rmpath('functions')

%% Display results
result_display(t,Lk,xk,xt,x_me,x_ul,x_vl)