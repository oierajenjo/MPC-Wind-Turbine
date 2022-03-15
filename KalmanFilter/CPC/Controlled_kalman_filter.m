clc
clear all
close all

%% Obtain all variables
variables_CPC
theta_f = 0;
[lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
u_b = [theta_f; K].*u_b;


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
% w_p = @(d) d(1)*pi/(2*W.L);
% ve = @(x,d) x(13) + d(1);
lamb = @(x,d) (x(1)*Ae.Rr)/(d(1));
% lamb = @(x,d) (x(1)*Ae.Rr-x(9))/(d(1)-x(7));
cp = @(x,d) cp_ct(lamb(x,d),x(10),cp_l,lambdaVec,pitchVec);
% cpy = @(x,d) cp_ct(lamb(x,d),x(10),cp_l,lambdaVec,pitchVec);
ct = @(x,d) cp_ct(lamb(x,d),x(10),ct_l,lambdaVec,pitchVec);

Tr = @(x,d) 0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^3*cp(x,d)/x(1);
Fx = @(x,d) 0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^2*ct(x,d);
Fy = @(x,d) (0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^3*cp(x,d)/x(1))/(2*Ae.Rr/3);

%% Drive train
f1 = @(x,d) (1-D.mu)*Tr(x,d)/(D.Jr+D.Jg) - x(12)/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3); % Tower foreafter velocity
f3 = @(x) -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.B*B.kx*x(6)/To.m + B.B*B.cx*x(7)/To.m; % Tower foreafter acceleration

f4 = @(x) x(5); % Tower edgewise velocity
f5 = @(x) 3*x(12)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.B*B.ky*x(8)/To.m + B.B*B.cy*x(9)/To.m; % Tower edgewise acceleration

%% Blades
f6 = @(x) x(7); % Blade foreafter velocity
f7 = @(x,d) Fx(x,d)/(B.B*B.m) + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(7)/B.m; % Blade foreafter acceleration

f8 = @(x) x(9); % Blade edgewise velocity
f9 = @(x,d) Fy(x,d)/(B.B*B.m) + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(8)/B.m - B.cy*x(9)/B.m; % Blade edgewise acceleration

%% Actuators
f10 = @(x) x(11); % Pitch velocity
f11 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(11) - Ac.omega^2*x(10); % Pitch acceleration
f12 = @(x,u) (u(2)*x(1)^2-x(12))/Ac.tau; % Torque change in time

%% Wind
% f13 = @(x,d) -w_p(d)*x(13); % Wind turbulence acceleration
% f14 = 0; % Mean wind acceleration


f = @(x,u,d) x + Ts*[f1(x,d); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x,d);...
    f8(x); f9(x,d); f10(x); f11(x,u); f12(x,u)]; % Nonlinear prediction

h = @(x,d) [x(1); f3(x); f5(x); B.l*B.m*f7(x,d); B.l*B.m*f9(x,d); ...
    D.eta*x(12)*x(1); d(1)];


% a = @(d) 1 - w_p(d)*Ts; % Euler
% a = @(d) exp(-w_p(d)*Ts); % Zero Order Hold
% sigma_t = @(d) W.ti*d(1)*sqrt((1-a(d)^2)/(1-a(d))^2);
% sigma_m = sqrt(Ts*W.q);
Q = @(d) diag([zeros(Lk-1,1); 0]); % Covariance matrix of the process noise

temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_root; M.sigma_root;...
    M.sigma_pow; M.sigma_vane].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Step 3: Initialize state and covariance
x = zeros(Lk, N); % Initialize size of state estimate for all k
% x(:,1) = [0]; % Set initial state estimate
x(:,1) = x_i;
P0 = 0.01*eye(Lk,Lk); % Set initial error covariance

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
yt = zeros(Yk, N); % Initialize size of output vector for all k

% Generate the true state values
for k = 2:N
    xt(:,k) = f(xt(:,k-1),u_b(:,k-1),d_b(:,k-1)) + Ts*n(d_b(:,k-1));
    yt(:,k-1) = h(xt(:,k-1),d_b(:,k-1)) + v(:,k-1);
end
for k=1:N
%     p(k) = vr(xt(:,k),d(:,k));
    frx(k) = Fx(xt(:,k),d_b(:,k))/(B.B*B.m);
    fry(k) = Fy(xt(:,k),d_b(:,k))/(B.B*B.m);
    tr(k) = (1-D.mu)*Tr(xt(:,k),d_b(:,k))/(D.Jr+D.Jg);
    cpl(k) = cp(xt(:,k),d_b(:,k));
    ctl(k) = ct(xt(:,k),d_b(:,k));
    la(k) = lamb(xt(:,k),d_b(:,k));
    %     pi(k) = vri(xt(:,k),1);
end

t = Ts*(1:N);

figure
plot(t,xt(1,:),t,d_b(3,:));
title("wr")
% xlim([1 50])
figure
plot(t,xt(2,:),t, -data.Data(:,224));
title("xt")
% xlim([1 50])
% figure
% plot(xt(3,:));
% title("xtdot")
figure
plot(t,xt(6,:));
title("xb")
% xlim([1 50])

figure
plot(t,xt(4,:),t, data.Data(:,225));
title("yt")
% xlim([1 50])
% figure
% plot(xt(5,:));
% title("ytdot")
figure
plot(t,xt(8,:));
title("yb")
% xlim([1 50])

figure
plot(t,d_b(1,:))
title("vr")
% xlim([1 50])
figure
plot(t,frx,t,fry)
title("Frx & Fry")


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
        chi_m(:,j) = f(chi_p(:,j),u_b(:,k-1),d_b(:,k-1));
    end
    
    x_m = chi_m*wm; % Calculate mean of predicted state
    % Calculate covariance of predicted state
    P_m = Q(d_b(:,k-1)); % A priori covariance estimate
    for i = 1:n_sigma_p
        P_m = P_m + wc(i)*(chi_m(:,i) - x_m)*(chi_m(:,i) - x_m)';
    end
    
    % Step 3: Observation Transformation
    % Propagate each sigma-point through observation
    % Initial velocity will be considered as 0, as we need it for
    % obtaining the acceleration
    psi_m = zeros(Yk,n_sigma_p);
    for j=1:n_sigma_p
        psi_m(:,j) = h(chi_m(:,j),d_b(:,k-1));
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
    xlabel('$Time (s)$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF');
    title(['$Error $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end

v_r = @(x,d) d(1,:)+x(13,:)-x(3,:);
figure(6)
plot(t,v_r(x,d),'b-',t,v_r(xt,d),'r-');
xlabel('Time [s]', 'FontSize', 14);
ylabel('Velocity [m/s]', 'FontSize', 14);
grid on;
legend('UKF', 'True');
title('Effective wind speed [v_r]', 'FontSize', 14);
set(gcf, 'PaperOrientation','landscape');
saveas(figure(6),'Figures/Kalman_ve.pdf');

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