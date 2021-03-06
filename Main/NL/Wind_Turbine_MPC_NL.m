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

[f,h,hz,Q,R] = system_IPC_NL(var,ct_l,cp_l,lambdaVec,pitchVec,Lk);

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
% xeq = x_i;
x_mpc = x_tv;
u_mpc = zeros(Uk, N);
% z_mpc = zeros(Zk, N);
% z_mpc(:,1) = z_i;

%% Reference trajectories
% ref_me = [Ac.omega_opt*ones(N+Hp,1), zeros(N+Hp,8), W.TSR*ones(N+Hp,3) , zeros(N+Hp,6), Ac.Pe_opt*ones(N+Hp,1)]';
% ref_me = Sz\ref_me;
ref_me = [Ac.omega_opt, zeros(1,8), W.TSR*ones(1,3) , zeros(1,6), Ac.Pe_opt];

%% MPC definition
MPCdefinition_NL

disp('Running Loop')
for k=1:N-1
    % Compute optimal control moves
    [u,nloptions] = nlmpcmove(nlobj,x_kf(:,k),uprev_mpc,ref_me,[],nloptions);
%     res = MPCobj({Sx\x_kf(:,k),uprev_mpc,ref_me(:,k+1:k+Hp)});

%     u_L = res{1};
%     u_temp = reshape(u_L, [Uk, length(u_L)/Uk]);
%     u_mpc(:,k) = u_temp(:,1);
%     uprev_mpc = u_temp(:,1);
    uprev_mpc = u;
    x_mpc(:,k+1) = nloptions.X0(1,:)';

%     x_L = res{2};
%     x_temp = Sx*reshape(x_L, [Lk, length(x_L)/Lk]);
%     x_mpc(:,k+1) = x_temp(:,1);
%     
%     z_L = res{3};
%     z_temp = reshape(z_L, [Zk, length(z_L)/Zk]);
%     z_mpc(:,k+1) = z_temp(:,1);
    
    %% Runge-Kutta 4th order method
    [x_tv(:,k+1),yt(:,k+1)] = RK4_NL(f,x_tv(:,k),uprev_mpc,h,n(x_tv(:,k)),v(:,k+1),Ts);

    %% Unscented Kalman Filter
    [x_kf(:,k+1),P,e(:,k+1)] = UKF_NL(f,h,Q,R,x_kf(:,k),yt(:,k+1),uprev_mpc,kal,P,Ts,v(:,k+1),n);
    
%     xeq = x_kf(:,k+1);
    if mod(k,30) == 0
         disp("Iteration: " + k);
    end
    toc
end
x_kf(end,:) = wrapToPi(x_kf(end,:))+pi;
x_tv(end,:) = wrapToPi(x_tv(end,:))+pi;

%% Display results
true_plots(Lk,yt,x_kf,x_ul,x_vl,t)
% result_display(t,Lk,x_kf,x_mpc,x_ul,x_vl)
save('working_MPC.mat')
rmpath('functions')