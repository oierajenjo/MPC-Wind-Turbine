clc
clear all
close all

%% Obtain all variables
Simple_Basic_variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

%% Controlled
% u_b = ones(1,N);
% 
% theta_f = 0;
% [lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
% K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
% u_b = [theta_f; K].*u_b;

%% Before filter execution
% Step 1: Define UT Scaling parameters and weight vectors
Lk = size(x_i,1); % Size of state vector
Yk = size(y_me,1); % Size of measured vector
Uk = size(u_b,1); % Size of imput vector

f1 = @(x) -3*x(3)/(2*To.H*To.m) - To.c*x(1)/To.m - To.k*x(2)/To.m; % Tower edgewise acceleration
f2 = @(x) x(1); % Tower edgewise velocity
f3 = @(x,u) (u(1)-x(3))/Ac.tau; % Torque change in time


f = @(x,u) [f1(x); f2(x); f3(x,u)]; % Nonlinear prediction
h = @(x) [f1(x)];

Q = diag(zeros(Lk,1)); % Covariance matrix of the process noise

temp = [M.sigma_acc].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Step 3: Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
% x(:,1) = [0]; % Set initial state estimate
xk(:,1) = x_i;

% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = sqrt(Q)*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
%w = zeros(1,N);
%v = zeros(1,N);
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k
yt(:,1) = y_me(:,1);

% Runge-Kutta 4th order method
for k = 2:N
    k_1 = f(xt(:,k-1),u_b(:,k-1));
    k_2 = f(xt(:,k-1)+0.5*Ts*k_1,u_b(:,k-1)+0.5*Ts);
    k_3 = f(xt(:,k-1)+0.5*Ts*k_2,u_b(:,k-1)+0.5*Ts);
    k_4 = f(xt(:,k-1)+Ts*k_3,u_b(:,k-1)+Ts);
    xt(:,k) = xt(:,k-1) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts;  % main equation
    
    yt(:,k) = h(xt(:,k)) + v(:,k);
end
% yt(:,N) = h(xt(:,N)) + v(:,N);

t = Ts*(1:N);
figure
plot(t,xt(2,:), t,data.Data(:,225));
title("yt")
legend(["Us" "Bladed"])
% xlim([1 50])
figure
plot(t,xt(1,:), t,data.Data(1,231));
title("ytdot")
legend(["Us" "Bladed"])

%% Kalman variables
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(Lk+kappa) - Lk;
n_sigma_p = 2*Lk + 1; % Number of sigma points
wm = ones(n_sigma_p,1)*1/(2*(Lk+lambda)); % Weight for transformed mean
wc = wm; % Weight for transformed covariance
wm(1) = lambda/(lambda+Lk);
wc(1) = lambda/(lambda+Lk) + 1 - alpha^2 + beta;

%% Execute Unscented Kalman Filter
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i;
% P0 = 0.01*eye(Lk,Lk); % Set initial error covariance
% P = P0; % Set first value of P to the initial P0
P0 = [0.01;0.01;0.01].^2;
P = diag(P0);
for k = 2:N
    % Step 1: Generate the sigma-points
    sP = chol(P,'lower'); % Calculate square root of error covariance
    % chi_p = "chi previous" = chi(k-1) % Untransformed sigma points
    chi_p = [xk(:,k-1), xk(:,k-1)*ones(1,Lk)+sqrt(Lk+lambda)*sP, ...
        xk(:,k-1)*ones(1,Lk)-sqrt(Lk+lambda)*sP]; % Untransformed sigma points
    
    % Step 2: Prediction Transformation
    % Propagate each sigma-point through prediction
    % chi_m = "chi minus" = chi(k|k-1)
    chi_m = zeros(Lk,n_sigma_p); % Transformed sigma points
    x_m = 0;
    for j=1:n_sigma_p
        chi_m(:,j) = chi_p(:,j) + Ts*f(chi_p(:,j),u_b(:,k-1));
        x_m = x_m + wm(j)*chi_m(:,j); % Calculate mean of predicted state
    end
%     x_m = chi_m*wm; % Calculate mean of predicted state
    % Calculate covariance of predicted state
    P_m = Q; % A priori covariance estimate
    for i = 1:n_sigma_p
        P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
    end
    
    % Step 3: Observation Transformation
    % Propagate each sigma-point through observation
    % Initial velocity will be considered as 0, as we need it for
    % obtaining the acceleration
    psi_m = zeros(Yk,n_sigma_p);
    y_m = 0;
    for j=1:n_sigma_p
        psi_m(:,j) = h(chi_m(:,j));
        y_m = y_m + wm(j)*psi_m(:,j); % Calculate mean of predicted output
    end
%     y_m = psi_m*wm; % Calculate mean of predicted output
    
    % Calculate covariance of predicted output
    % and cross-covariance between state and output
    Pyy = R;
    Pxy = 0;
    for i = 1:n_sigma_p
        Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
        Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
    end
    
    % Step 4: Measurement Update
    K = Pxy/Pyy; % Calculate Kalman gain
    xk(:,k+1) = x_m + K*(yt(:,k+1) - y_m); % Update state estimate
    P = P_m - K*Pyy*K'; % Update covariance estimate
end

%% Display results
t = 1:N;
for i = 1:Lk
    figure
    subplot(1,2,1);
    %     plot(t,xt(i,:),'r-','LineWidth', 2);
    plot(t,xk(i,:),'b-', t,xt(i,:),'r-');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF','True');
    title(['$Estimations $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
    
    subplot(1,2,2);
    plot(t,xk(i,:)-xt(i,:),'b-');
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF');
    title(['$Error $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end   
