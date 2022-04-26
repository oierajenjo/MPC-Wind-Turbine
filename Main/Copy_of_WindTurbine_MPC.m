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


[f,h,Q,R] = system_IPC(var,ct_l,cp_l,lambdaVec,pitchVec,Lk);

% Step 3: Initialize state and covariance
% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(x) sqrt(Q(x))*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
% v = zeros(Yk, N);

%% Initialization
% Initialize matrices
x_tv = zeros(Lk, N); % Initialize size of true state for all k
x_tv(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k

% Initialize state and covariance
x_kf = zeros(Lk, N); % Initialize size of state estimate for all k
x_kf(:,1) = x_i;
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
xmpc = zeros(Lk, N); % Initialize size of true state for all k
xmpc(:,1) = x_i; % Set true initial state


disp('Running Loop')

%% Reference trajectories
ref_me = [Ac.omega_opt*ones(N,1) zeros(N,14) Ac.Pe_opt*ones(N,1)]';

for k=1:N-1
    %% Linearized system
    A1 = [zeros(1,11), -a1*ones(1,3);
        zeros(1,2), 1, zeros(1,11);
        0, -b3, -b4, zeros(1,2), b1*ones(1,3), b2*ones(1,3), zeros(1,3);
        zeros(1,4), 1, zeros(1,9);
        zeros(1,3), -c4, -c5, zeros(1,6), c2*ones(1,3);
        zeros(3,8), eye(3), zeros(3,3);
        d10(xeq,0), d3, d4(xeq,0), zeros(1,2), -d1, zeros(1,2), d2(xeq,0), zeros(1,5);
        d10(xeq,1), d3, d4(xeq,1), zeros(1,3), -d1, zeros(1,2), d2(xeq,1), zeros(1,4);
        d10(xeq,2), d3, d4(xeq,2), zeros(1,4), -d1, zeros(1,2), d2(xeq,2), zeros(1,3);
        zeros(3,14)];
    
    A2 = [zeros(1,9), -a2, zeros(1,3);
        zeros(3,13);
        c3*ones(1,3), zeros(1,6), -c1, zeros(1,3);
        zeros(3,13);
        d8(xeq,0), zeros(1,2), d7(xeq,0), zeros(1,6), d6(xeq,0), d5(xeq,0), d9(xeq,0);
        0, d8(xeq,1), zeros(1,2), d7(xeq,1), zeros(1,5), d6(xeq,1), d5(xeq,1), d9(xeq,1);
        0, 0, d8(xeq,2), zeros(1,2), d7(xeq,2), zeros(1,4), d6(xeq,2), d5(xeq,2), d9(xeq,2);
        eye(3), zeros(3,10)];
    
    A3 = [-e9(xeq,0), 0, -e7(xeq,0), e3, e4, zeros(1,3), -e8(xeq,0), 0, 0, -e1, 0, 0;
        -e9(xeq,1), 0, -e7(xeq,1), e3, e4, zeros(1,4), -e8(xeq,1), 0, 0, -e1, 0;
        -e9(xeq,2), 0, -e7(xeq,2), e3, e4, zeros(1,5), -e8(xeq,2), 0, 0, -e1;
        zeros(6,14);
        g2(xeq,u_b(:,k)), zeros(1,13);
        zeros(2,14);
        1, zeros(1,13)];
    
    A4 = [-e2(xeq,0), zeros(1,2), -e11(xeq,0), zeros(1,6), -e6(xeq,0), -e5(xeq,0), -e10(xeq,0);
        0, -e2(xeq,1), zeros(1,2), -e11(xeq,1), zeros(1,5), -e6(xeq,1), -e5(xeq,1), -e10(xeq,1);
        0, 0, -e2(xeq,2), zeros(1,2), -e11(xeq,2), zeros(1,4), -e6(xeq,2), -e5(xeq,2), -e10(xeq,2);
        zeros(3,6), eye(3), zeros(3,4);
        zeros(3), -f1*eye(3), -f2*eye(3), zeros(3,4);
        zeros(1,9), -g1, zeros(1,3);
        zeros(1,10), -W.w_p(xeq), zeros(1,2);
        zeros(2,13)];
    
    Ampc = [A1 A2;
        A3 A4];
    % Ampc = eye(Lk); % State Matrix
    
    Bmpc = [zeros(3,Lk-7), f1*eye(3), zeros(3,4);
        zeros(1,Lk)]';
%     Bmpc = [zeros(3,Lk-7), f1*eye(3), zeros(3,4);
%         zeros(1,Lk-4), g1, zeros(1,3)]';
    % Bmpc = Ts*eye(Lk,Uk); % Input Matrix
    
    % CHANGE
    Cmpc = [1, zeros(1,Lk-1);
        0, 1, zeros(1,Lk-2);
        zeros(1,4), 1, zeros(1,Lk-5);
        zeros(3,8), eye(3), zeros(3,Lk-8-3);
        zeros(3,14), eye(3), zeros(3,Lk-14-3);
        zeros(3,Lk-7-3), eye(3), zeros(3,7);
        zeros(3,Lk-4-3), eye(3), zeros(3,4);
        q2(xeq), zeros(1,Lk-5), q1(xeq), zeros(1,3)];
%     xmpc(:,k+1) = xmpc(:,k) + Ts*(Ampc*xmpc(:,k) + Bmpc*u_b(:,k));
%     ympc(:,k+1) = Cmpc*xmpc(:,k+1);

    k_1 = Ampc*xmpc(:,k) + Bmpc*u_b(:,k);
    k_2 = Ampc*(xmpc(:,k)+0.5*Ts*k_1) + Bmpc*(u_b(:,k)+0.5*Ts);
    k_3 = Ampc*(xmpc(:,k)+0.5*Ts*k_2) + Bmpc*(u_b(:,k)+0.5*Ts);
    k_4 = Ampc*(xmpc(:,k)+Ts*k_3) + Bmpc*(u_b(:,k)+Ts);
    xmpc(:,k+1) = xmpc(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts; 
    ympc(:,k+1) = Cmpc*xmpc(:,k+1);
    
    %% Runge-Kutta 4th order method
    % disp('Running True Values')
    [x_tv(:,k+1),yt(:,k+1)] = RK4(f,x_tv(:,k),u_b(:,k),h,n(x_tv(:,k)),v(:,k+1),Ts);
    
    xeq = x_tv(:,k+1);
end
xmpc(end,:) = wrapToPi(xmpc(end,:))+pi;

%% Display results
figure
plot(t,xmpc(1,:),t,x_tv(1,:));
legend('$\omega_r$ (linear)','$\omega_r$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Rotor speed');
ylabel('$\omega_r$ [rad/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(3,:),t,x_tv(3,:));
legend('$\dot{x}_{t}$ (linear)','$\dot{x}_{t}$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower fore-aft velocity');
ylabel('$\dot{x}_{t}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(5,:),t,x_tv(5,:));
legend('$\dot{y}_{t}$ (linear)','$\dot{y}_{t}$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower sidewards velocity');
ylabel('$\dot{y}_{t}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

% figure
% plot(t,xmpc(4,:)',t,x_tv(4,:));
% legend('$M_{y_1}$ (sim.)','$M_{y_1}$ (Bladed)','Interpreter','latex')
% % ylim([-2.6 2.6]);
% title('Flapwise root bending moment (blade 1)');
% ylabel('$M_{y_1}$ [Nm]','Interpreter','latex')
% xlabel('time [s]','Interpreter','latex')

% figure
% plot(t,yt(7,:)',t,y_me(7,:));
% legend('$M_{x_1}$ (sim.)','$M_{x_1}$ (Bladed)','Interpreter','latex')
% % ylim([-2.6 2.6]);
% title('Edgewise root bending moment (blade 1)');
% ylabel('$M_{x_1}$ [Nm]','Interpreter','latex')
% xlabel('time [s]','Interpreter','latex')

figure
plot(t,ympc(16,:),t,yt(13,:));
legend('$P_e$ (linear)','$P_e$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Electrical power');
ylabel('$P_e$ [W]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

% figure
% plot(t,xmpc(11,:)',t,y_me(11,:));
% legend('$v_r$ (sim.)','$v_r$ (Bladed)','Interpreter','latex')
% % ylim([-2.6 2.6]);
% title('Relative wind speed at the rotor');
% ylabel('$v_r$ [m/s]','Interpreter','latex')
% xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(2,:),t,x_tv(2,:));
legend('$x_t$ (linear)','$x_t$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower fore-aft deflection');
ylabel('$x_t$ [m]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(4,:),t,x_tv(4,:));
legend('$y_t$ (linear)','$y_t$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower sidewards deflection');
ylabel('$y_t$ [m]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(6,:),t,x_tv(6,:));
legend('$x_{b_1}$ (linear)','$x_{b_1}$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Flapwise deflection (blade 1)');
ylabel('$x_{b_1}$ [m]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(9,:),t,x_tv(9,:));
legend('$\dot{x}_{b_1}$ (linear)','$\dot{x}_{b_1}$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Flapwise velocity (blade 1)');
ylabel('$\dot{x}_{b_1}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(12,:),t,x_tv(12,:));
legend('$y_{b_1}$ (linear)','$y_{b_1}$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Edgewise deflection (blade 1)');
ylabel('$y_{b_1}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(15,:),t,x_tv(15,:));
legend('$\dot{y}_{b_1}$ (linear)','$\dot{y}_{b_1}$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Edgewise velocity (blade 1)');
ylabel('$\dot{y}_{b_1}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(18,:),t,x_tv(18,:));
legend('$\theta_1 (linear)$','$\theta_1$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Pitch angle (blade 1)');
ylabel('$\theta_1$ [rad]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(22,:),t,x_tv(22,:));
legend('$\dot{\theta}_1$ (linear)','$\dot{\theta}_1$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Pitch rate (blade 1)');
ylabel('$\dot{\theta}_1$ [rad/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(24,:),t,x_tv(24,:));
legend('$T_g$ (linear)','$T_g$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Generator torque');
ylabel('$T_g$ [Nm]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(25,:),t,x_tv(25,:));
legend('$v_t$ (linear)','$v_t$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Turbulent wind speed');
ylabel('$v_t$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xmpc(26,:),t,x_tv(26,:));
legend('$v_m$ (linear)','$v_m$ (non-linear)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Mean wind speed');
ylabel('$v_m$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')


% true_plots(ympc,y_me,xmpc,data,t)
% result_display(t,Lk,x_kf,x_tv,x_me,x_ul,x_vl)
rmpath('functions')