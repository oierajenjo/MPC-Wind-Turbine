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
% kAns = 'No';
% kAns = convertCharsToStrings(kAns);
% if kAns=="Cancel"
%     return
% elseif kAns == "Yes"
%     disp("Measured values")
% elseif kAns == "No"
%     disp("True values")
% end

% u_b = ones(4,N);
% theta_f = 0;
% [lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
% K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
% u_b = [theta_f; theta_f; theta_f; K].*u_b;

[f,h,Q,R] = system_IPC(var,ct_l,cp_l,lambdaVec,pitchVec,Lk);

% Step 3: Initialize state and covariance
% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(x) sqrt(Q)*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
% v = zeros(Yk, N);

%% Initialization
% Initialize matrices
x_tv = zeros(Lk, N); % Initialize size of true state for all k
x_tv(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k

% Initialize state and covariance
x_kf = x_tv; % Initialize size of state estimate for all k
P0 = [M.sigma_enc; M.sigma_tdef; M.sigma_tvel; M.sigma_tdef; M.sigma_tvel;...
    M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel; M.sigma_bvel;...
    M.sigma_bvel; M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel;...
    M.sigma_bvel; M.sigma_bvel; M.sigma_pit; M.sigma_pit; M.sigma_pit;...
    M.sigma_pitvel; M.sigma_pitvel; M.sigma_pitvel; M.sigma_pow;...
    M.sigma_vane; 0.01; M.sigma_azim].^2;
P0 = diag(P0);
P = P0;
e = zeros(Yk, N);

uprev_mpc = u_b(:,1);
xeq = x_i;
x_mpc = x_tv;
u_mpc = zeros(Uk, N);

%% Reference trajectories
ref_me = [Ac.omega_opt*ones(N+Hp,1), zeros(N+Hp,8), W.TSR*ones(N+Hp,3) , zeros(N+Hp,6), Ac.Pe_opt*ones(N+Hp,1)]';
ref_me = Sz\ref_me;

disp('Running Loop')
for k=1:N-1
    %% MPC
    MPCdefinition
    res = MPCobj({Sx\x_kf(:,k),uprev_mpc,ref_me(:,k+1:k+Hp)});

    u_L = res{1};
    u_temp = reshape(u_L, [Uk, length(u_L)/Uk]);
    u_mpc(:,k) = u_temp(:,1);
    uprev_mpc = u_temp(:,1);

    x_L = res{2};
    x_temp = Sx*reshape(x_L, [Lk, length(x_L)/Lk]);
    x_mpc(:,k+1) = x_temp(:,1);
    
    %% Runge-Kutta 4th order method
%     disp('Running True Values')
    [x_tv(:,k+1),yt(:,k+1)] = RK4(f,x_tv(:,k),uprev_mpc,h,n(x_tv(:,k)),v(:,k+1),Ts);

    %% Unscented Kalman Filter
%     if kAns == "Yes"
%         disp('Running Kalman Filter')
%         [x_kf(:,k+1),P,e(:,k+1)] = UKF(f,h,Q,R,x_kf(:,k),y_me(:,k+1),uprev_mpc,kal,P,Ts,v(:,k+1),n);
%     elseif kAns == "No"
%         disp('Running Kalman Filter')
    [x_kf(:,k+1),P,e(:,k+1)] = UKF(f,h,Q,R,x_kf(:,k),yt(:,k+1),uprev_mpc,kal,P,Ts,v(:,k+1),n);
%     end
    
    xeq = x_tv(:,k+1);
    if mod(k,30) == 0
         disp("Iteration: " + k);
    end
end
x_kf(end,:) = wrapToPi(x_kf(end,:))+pi;
x_tv(end,:) = wrapToPi(x_tv(end,:))+pi;

%% Display results
true_plots(yt,y_me,x_tv,data,t)
result_display(t,Lk,x_kf,x_tv,x_ul,x_vl)
rmpath('functions')