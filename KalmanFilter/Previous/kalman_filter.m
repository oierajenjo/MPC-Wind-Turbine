clc
clear all
close all

%% Obtain all variables
run("variables.m")

%% Before filter execution
% System properties
N = data.Channels.Scans; % Number of time steps for filter
%N1 = 20; % Station 1 North coordinate
%E1 = 0; % Station 1 East coordinate
%N2 = 0; % Station 2 North coordinate
%E2 = 20; % Station 2 East coordinate

% Step 1: Define UT Scaling parameters and weight vectors
Lk = 9; % Size of state vector
Yk = 5; % Size of measured vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(Lk+kappa) - Lk;
n_sigma_p = 2*Lk + 1; % Number of sigma points
wm = ones(n_sigma_p,1)*1/(2*(Lk+lambda)); % Weight for transformed mean
wc = wm; % Weight for transformed covariance
wm(1) = lambda/(lambda+Lk);
wc(1) = lambda/(lambda+Lk) + 1 - alpha^2 + beta;

% Step 2: Define noise assumptions
wp = @(d) d(1)*pi/(2*L);
ve = @(x,d) x(9)+d(1);
Tr = @(x,d) 0.5*rho*(ve(x,d)-x(3))^3*Ar*cp_ct(x(1)*r/(ve(x,d)-x(3)),x(6),cp_l,lambdaVec,pitchVec)/x(1);
Fr = @(x,d) 0.5*rho*(ve(x,d)-x(3))^2*Ar*cp_ct(x(1)*r/(ve(x,d)-x(3)),x(6),ct_l,lambdaVec,pitchVec);

%% Drive train
f1 = @(x,d) ((1-mu_d)*Tr(x,d) - x(8))/(Jr+Jg);

%% Tower
f2 = @(x) x(3);
f3 = @(x,d) (Fr(x,d)-ct*x(3)-kt*x(2))/mt;
f4 = @(x) x(5);
f5 = @(x) ((3/(2*H))*x(8)-ct*x(5)-kt*x(4))/mt;

%% Actuators
f6 = @(x) x(7); % Pitch velocity
f7 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(7) - Ac.omega^2*x(6); % Pitch acceleration
f8 = @(x,u) (u(2)-x(8))/Ac.tau; % Torque change in time

%% Wind
f9 = @(x,d) -wp(d)*x(9);
% f10 = 0;

f = @(x,u,d) x + Ts*[f1(x,d); f2(x); f3(x,d); f4(x); f5(x); f6(x); f7(x,u); f8(x,u); f9(x,d)]; % Nonlinear prediction

% h = @(x) (x);
h = @(x,d) [x(1); f3(x,d); f5(x); ve(x,d)-x(3); eta_g*x(8)*x(1)];

a = @(d) 1- wp(d)*Ts; % Euler
% a = @(x) exp(-(x(5)*pi/(2*L))*Ts); % Zero Order Hold
sigma_t = @(d) ti*d(1)*sqrt((1-a(d)^2)/(1-a(d))^2);
sigma_m = sqrt(Ts*q);
Q = @(d) diag([zeros(Lk-1,1); sigma_t(d)^2*wp(d)^2]); % Covariance matrix of the process noise
temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_pow; M.sigma_vane].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Step 3: Initialize state and covariance
x = zeros(Lk, N); % Initialize size of state estimate for all k
% x(:,1) = [0]; % Set initial state estimate
x(:,1) = x_i;
P0 = eye(Lk,Lk); % Set initial error covariance

% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(d) sqrt(Q(d))*randn(Lk, 1); % Generate random process noise (from assumed Q) 
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
%w = zeros(1,N);
%v = zeros(1,N);
xt = zeros(Lk, N); % Initialize size of true state for all k
% xt(:,1) = zeros(Lk,1) + sqrt(P0)*randn(Lk,1); % Set true initial state
xt(:,1) = x_i; % Set true initial state
y = zeros(Yk, N); % Initialize size of output vector for all k
% y(:,1) = h(xt(:,1),u(:,1)) + v(:,1);
% Generate the true state values
for k = 2:N
    xt(:,k) = f(xt(:,k-1),u(:,k-1),d(:,k-1)) + Ts*n(d(:,k-1));
    y(:,k-1) = h(xt(:,k-1),d(:,k-1)) + v(:,k-1);
end
for k=1:N
    p(k) = ve(xt(:,k),d(:,k));
    fr(k) = Fr(xt(:,k),d(:,k));
    tr(k) = (1-mu_d)*Tr(xt(:,k),d(:,k));
%     pi(k) = vri(xt(:,k),1);
end
figure
plot(xt(1,:));
title("wr")
figure
plot(xt(2,:));
title("xt")
figure
plot(xt(3,:));
title("xtdot")
figure
plot(xt(4,:));
title("yt")
figure
plot(xt(5,:));
title("ytdot")
figure
plot(p(:))
title("ve")
figure
plot(1:N,xt(8,:),1:N,tr(:))
title("Tg & Tr")

figure
plot(3*xt(8,:)/(2*H*mt))
title("Tg Applied")

figure
plot(fr(:))
title("Fr")
% %% Initialize and run EKF for comparison
% xe = zeros(Lk,N);
% xe(:,1) = x(:,1);
% P = P0;
% for k = 2:N
%     % Prediction
%     x_m = f(xe(:,k-1),u(:,k-1));
%     F = 1-Ts^2*sin(Ts*k);
%     % F = -T*sin(T*k);
%     P_m = F*P*F' + Q;
% 
%     % Observation
%     y_m = x_m(1);
% 
%     H = 1;
%     %H = -T*(sin(T*x_m(1)) + T*cos(T*x_m(1)));
% 
%     % Measurement Update
%     K = P_m*H'/(H*P_m*H' + R); % Calculate Kalman gain
%     xe(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
%     P = (eye(Lk)-K*H)*P_m; % Update covariance estimate
% end

%% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0
for k = 2:N
    % Step 1: Generate the sigma-points
    sP = chol(P,'lower'); % Calculate square root of error covariance
    % chi_p = "chi previous" = chi(k-1) % Untransformed sigma points
    chi_p = [x(:,k-1), x(:,k-1)*ones(1,Lk)+sqrt(Lk+lambda)*sP, ...
             x(:,k-1)*ones(1,Lk)-sqrt(Lk+lambda)*sP]; % Untransformed sigma points

    % Step 2: Prediction Transformation
    % Propagate each sigma-point through prediction
    % chi_m = "chi minus" = chi(k|k-1)
    chi_m = zeros(Lk,n_sigma_p); % Transformed sigma points
    for j=1:n_sigma_p
        chi_m(:,j) = f(chi_p(:,j),u(:,k-1),d(:,k-1));
    end
    
    x_m = chi_m*wm; % Calculate mean of predicted state
    % Calculate covariance of predicted state
    P_m = Q(d(:,k-1)); % A priori covariance estimate
    for i = 1:n_sigma_p
        P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
    end

    % Step 3: Observation Transformation
    % Propagate each sigma-point through observation
    % Initial velocity will be considered as 0, as we need it for
    % obtaining the acceleration
    psi_m = zeros(Yk,n_sigma_p);
    for j=1:n_sigma_p
        psi_m(:,j) = h(chi_m(:,j),d(:,k-1));
    end
    y_m = psi_m*wm; % Calculate mean of predicted output
    
    % Calculate covariance of predicted output
    % and cross-covariance between state and output
    Pyy = R;
    %Pxy = zeros(L,2);
    Pxy = 0;
    for i = 1:n_sigma_p
        Pyy = Pyy + wc(i)*(psi_m(:,i) - y_m)*(psi_m(:,i) - y_m)';
        Pxy = Pxy + wc(i)*(chi_m(:,i) - x_m)*(psi_m(:,i) - y_m)';
    end

    % Step 4: Measurement Update
    K = Pxy/Pyy; % Calculate Kalman gain
    x(:,k) = x_m + K*(y(:,k) - y_m); % Update state estimate
    P = P_m - K*Pyy*K'; % Update covariance estimate
end

%% Display results
t = Ts*(1:N);
for i = 1:Lk
    figure(i)
    subplot(1,2,1); 
%     plot(t,xt(i,:),'r-','LineWidth', 2);
    plot(t,x(i,:),'b-', t,xt(i,:),'r-');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14); 
    grid on; 
    legend('UKF','True');
    title(['$Estimations $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
    
    subplot(1,2,2); 
    plot(t,x(i,:)-xt(i,:),'b-');
    xlabel('$Time (s)$', 'Interpreter', 'latex', 'FontSize', 14); 
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14); 
    grid on;
    legend('UKF');
    title(['$Error $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end

v_r = @(x,d) x(4,:)+d(1,:)-x(2,:);
figure(6)
plot(t,v_r(x,d),'b-',t,v_r(xt,d),'r-');
xlabel('Time [s]', 'FontSize', 14);
ylabel('Velocity [m/s]', 'FontSize', 14);
grid on;
legend('UKF', 'True');
title('Effective wind speed [v_r]', 'FontSize', 14);
set(gcf, 'PaperOrientation','landscape');
% saveas(figure(6),'Figures/Kalman_ve.pdf');

function res = cp_ct(la,be,cl,lambdaVec,pitchVec)
    [~,i_la] = min(abs(lambdaVec-abs(la)));
    [~,i_be] = min(abs(pitchVec-be));
    res = cl(i_la,i_be);
end