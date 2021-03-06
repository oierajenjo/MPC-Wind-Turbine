clc
clear all
close all

%% Obtain all variables
Basic_variables_CPC
% u_b = ones(2,N);
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
w_p = @(x) x(10)*pi/(2*W.L);
ve = @(x) x(10) + x(9);
vr = @(x) ve(x) - x(3);

lamb = @(x) (x(1)*Ae.Rr)/vr(x);
cp = @(x) cp_ct(lamb(x),x(6),cp_l,lambdaVec,pitchVec);
ct = @(x) cp_ct(lamb(x),x(6),ct_l,lambdaVec,pitchVec);

Tr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x))^3*cp(x)/x(1);
Fx = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x))^2*ct(x);
Fy = @(x) (0.5*Ae.rho*Ae.Ar*(vr(x))^3*cp(x)*3)/(2*x(1)*Ae.Rr);

%% Drive train
f1 = @(x) (1-D.mu)*Tr(x)/(D.Jr+D.Jg) - x(8)/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3); % Tower foreafter velocity
f3 = @(x) Fx(x)/To.m - To.c*x(3)/To.m - To.k*x(2)/To.m; % Tower foreafter acceleration

f4 = @(x) x(5); % Tower edgewise velocity
f5 = @(x) -3*x(8)/(2*To.H*To.m) - To.c*x(5)/To.m - To.k*x(4)/To.m; % Tower edgewise acceleration

%% Actuators BIEN
f6 = @(x) x(7); % Pitch velocity
f7 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(7) - Ac.omega^2*x(6); % Pitch acceleration
f8 = @(x,u) (u(2)-x(8))/Ac.tau; % Torque change in time

%% Wind
f9 = @(x) -w_p(x)*x(9); % Wind turbulence acceleration
f10 = 0; % Mean wind acceleration


f = @(x,u) [f1(x); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x,u);...
    f8(x,u); f9(x); f10]; % Nonlinear prediction

h = @(x) [x(1); f3(x); f5(x); D.eta*x(8)*x(1); vr(x)];


a = @(x) 1 - w_p(x)*Ts; % Euler
% a = @(x) exp(-w_p(x)*Ts); % Zero Order Hold
sigma_t = @(x) W.ti*x(10)*sqrt((1-a(x)^2)/(1-a(x))^2);
sigma_m = sqrt(Ts*W.q);
Q = @(x) diag([zeros(Lk-2,1); sigma_t(x)*w_p(x)^2; sigma_m]); % Covariance matrix of the process noise

temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_pow; M.sigma_vane].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(x) sqrt(Q(x))*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k

% % Generate the true state values
% for k = 2:N
%     xt(:,k) = xt(:,k-1) + Ts*(f(xt(:,k),u_b(:,k))+f(xt(:,k-1),u_b(:,k-1)))/2 ...
%     + Ts*n(xt(:,k-1));
%     yt(:,k-1) = h(xt(:,k-1)) + v(:,k-1);
% end

% Runge-Kutta 4th order method
for k = 1:N-1  
    k_1 = f(xt(:,k),u_b(:,k));
    k_2 = f(xt(:,k)+0.5*Ts*k_1,u_b(:,k)+0.5*Ts);
    k_3 = f(xt(:,k)+0.5*Ts*k_2,u_b(:,k)+0.5*Ts);
    k_4 = f(xt(:,k)+Ts*k_3,u_b(:,k)+Ts);
    xt(:,k+1) = xt(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n(xt(:,k));  % main equation
    
    yt(:,k) = h(xt(:,k)) + v(:,k);
end
yt(:,N) = h(xt(:,N)) + v(:,N);

for k=1:N
    fx(k) = Fx(xt(:,k))/(To.m);
    fy(k) = Fy(xt(:,k))/(To.m);
    tr(k) = (1-D.mu)*Tr(xt(:,k))/(D.Jr+D.Jg);
    p(k) = ve(xt(:,k));
end

t = Ts*(1:N);
figure
plot(t,xt(1,:), t,data.Data(:,10));
title("wr")
legend(["Us" "Bladed"])
% xlim([1 50])
figure
plot(t,xt(2,:), t,-data.Data(:,224));
title("xt")
legend(["Us" "Bladed"])
% xlim([1 50])
% figure
% plot(xt(3,:));
% title("xtdot")

figure
plot(t,xt(4,:), t,data.Data(:,225));
title("yt")
legend(["Us" "Bladed"])
% xlim([1 50])
% figure
% plot(xt(5,:));
% title("ytdot")

figure
plot(t,yt(5,:), t,data.Data(:,54))
title("vr")
legend(["Us" "Bladed"])
% xlim([1 50])

figure
plot(t,fx, t,fy)
title("Fx & Fy")
% legend(["Us" "Bladed"])
% xlim([1 50])

% Execute Unscented Kalman Filter
% Step 3: Initialize state and covariance
x = zeros(Lk, N); % Initialize size of state estimate for all k
x(:,1) = x_i;
P0 = 0.01*eye(Lk,Lk); % Set initial error covariance
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
        chi_m(:,j) = f(chi_p(:,j),u_b(:,k-1));
    end
    
    x_m = chi_m*wm; % Calculate mean of predicted state
    % Calculate covariance of predicted state
    P_m = Q(x(:,k-1)); % A priori covariance estimate
    for i = 1:n_sigma_p
        P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
    end
    
    % Step 3: Observation Transformation
    % Propagate each sigma-point through observation
    % Initial velocity will be considered as 0, as we need it for
    % obtaining the acceleration
    psi_m = zeros(Yk,n_sigma_p);
    for j=1:n_sigma_p
        psi_m(:,j) = h(chi_m(:,j));
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
    x(:,k) = x_m + K*(yt(:,k) - y_m); % Update state estimate
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
  ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF');
    title(['$Error $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end

v_r = @(x) x(14,:)+x(13,:)-x(3,:);
figure(6)
plot(t,v_r(x),'b-',t,v_r(xt,d),'r-');
xlabel('Time [s]', 'FontSize', 14);
ylabel('Velocity [m/s]', 'FontSize', 14);
grid on;
legend('UKF', 'True');
title('Effective wind speed [v_r]', 'FontSize', 14);
set(gcf, 'PaperOrientation','landscape');
saveas(figure(6),'Figures/Kalman_ve.pdf');   xlabel('$Time (s)$', 'Interpreter', 'latex', 'FontSize', 14);
   

function [la,res] = cp_max(be,cl,lambdaVec,pitchVec)
[~,i_be] = min(abs(pitchVec-be));
l_c = cl(:,i_be);
[res,i_la] = max(l_c);
la = lambdaVec(i_la);
end

function res = cp_ct(la,be,cl,lambdaVec,pitchVec)
[~,i_la] = min(abs(lambdaVec-abs(la)));
[~,i_be] = min(abs(pitchVec-be));
res = cl(i_la,i_be);
end