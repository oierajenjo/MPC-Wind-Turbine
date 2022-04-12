clc
clear all
close all
rng(1);

%% Obtain all variables
variables_IPC
load('BladedFiles\performancemap_data.mat')
Constant_variables
MPCconstants
addpath('functions');

% UNCOMMENT
% kAns = questdlg('Execute Kalman Filter with Bladed measurements', ...
% 	'Kalman Filter');
kAns = 'No';
kAns = convertCharsToStrings(kAns);
if kAns=="Cancel"
    return
elseif kAns == "Yes"
    disp("Measured values")
elseif kAns == "No"
    disp("True values")
end

u_b = ones(4,N);
theta_f = 0;
[lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
u_b = [theta_f; theta_f; theta_f; K].*u_b;


[f,h,Q,R] = system_IPC(var,Ts,ct_l,cp_l,lambdaVec,pitchVec,Lk);

% Step 3: Initialize state and covariance
% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(x) sqrt(Q(x))*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
% v = zeros(Yk, N);

%% Initialization
% Initialize matrices
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k

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
P = P0;
e = zeros(Yk, N);

disp('Running Loop')
for k=1:N-1
    %% MPC
    

    %% Runge-Kutta 4th order method
%     disp('Running True Values')
    [xt(:,k+1),yt(:,k+1)] = RK4(f,xt(:,k),u_b(:,k),h,n(xt(:,k)),v(:,k+1),Ts);

    %% Unscented Kalman Filter
    if kAns == "Yes"
%         disp('Running Kalman Filter')
        [xk(:,k+1),P,e(:,k+1)] = UKF(f,h,Q,R,xk(:,k),y_me(:,k+1),u_b(:,k),kal,P,Ts,v(:,k),n);
    elseif kAns == "No"
%         disp('Running Kalman Filter')
        [xk(:,k+1),P,e(:,k+1)] = UKF(f,h,Q,R,xk(:,k),yt(:,k+1),u_b(:,k),kal,P,Ts,v(:,k),n);
    end
end
xk(end,:) = wrapToPi(xk(end,:))+pi;
xt(end,:) = wrapToPi(xt(end,:))+pi;

%% Display results
true_plots(yt,y_me,xt,data,t)
result_display(t,Lk,xk,xt,x_me,x_ul,x_vl)
rmpath('functions')