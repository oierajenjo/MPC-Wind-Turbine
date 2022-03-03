clc
clear all
close all

%% Obtain all variables
variables_CPC

%% Before filter execution
% System properties
%N1 = 20; % Station 1 North coordinate
%E1 = 0; % Station 1 East coordinate
%N2 = 0; % Station 2 North coordinate
%E2 = 20; % Station 2 East coordinate

% Step 1: Define UT Scaling parameters and weight vectors
Lk = size(x_i,1); % Size of state vector
Yk = size(y_me,1); % Size of measured vector
Uk = size(u,1); % Size of imput vector
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
w_p = @(d) d(1)*pi/(2*W.L);
ve = @(x,d) x(13) + d(1);
Tr = @(x,d) 0.5*Ae.rho*Ae.Ar*(ve(x,d)-x(3))^3*cp_ct(x(1)*Ae.Rr/(ve(x,d)-x(3)),x(10),cp_l,lambdaVec,pitchVec)/x(1);
Fr = @(x,d) 0.5*Ae.rho*Ae.Ar*(ve(x,d)-x(3)-x(7))^2*cp_ct(x(1)*Ae.Rr/(ve(x,d)-x(3)),x(10),ct_l,lambdaVec,pitchVec);

%% Drive train
f1 = @(x,d) (1-D.mu)*Tr(x,d)/(D.Jr+D.Jg) - x(12)/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3); % Tower foreafter velocity
f3 = @(x) -(B.B*B.kx + T.k)*x(2)/T.m - (B.B*B.cx + T.c)*x(3)/T.m + B.B*B.kx*x(6)/T.m + B.B*B.cx*x(7)/T.m; % Tower foreafter acceleration

f4 = @(x) x(5); % Tower edgewise velocity
f5 = @(x) 3*x(12)/(2*T.H*T.m) - (B.B*B.ky + T.k)*x(4)/T.m - (B.B*B.cy + T.c)*x(5)/T.m + B.B*B.ky*x(8)/T.m + B.B*B.cy*x(9)/T.m; % Tower edgewise acceleration

%% Blades
f6 = @(x) x(7); % Blade foreafter velocity
f7 = @(x,d) Fr(x,d)/(B.B*B.m) + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(7)/B.m; % Blade foreafter acceleration

f8 = @(x) x(9); % Blade edgewise velocity
f9 = @(x) B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(8)/B.m - B.cy*x(9)/B.m; % Blade edgewise acceleration

%% Actuators BIEN
f10 = @(x) x(11); % Pitch velocity
f11 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(11) - Ac.omega^2*x(10); % Pitch acceleration
f12 = @(x,u) (u(2)-x(12))/Ac.tau; % Torque change in time

%% Wind
f13 = @(x,d) -w_p(d)*x(13); % Wind turbulence acceleration
% f14 = 0; % Mean wind acceleration


f = @(x,u,d) x + Ts*[f1(x,d); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x,d);...
    f8(x); f9(x); f10(x); f11(x,u); f12(x,u); f13(x,d)]; % Nonlinear prediction

h = @(x,d) [x(1); f3(x); f5(x); B.l*B.m*f7(x,d); B.l*B.m*f9(x); ...
    D.eta*x(12)*x(1); ve(x,d)-x(3)];


a = @(d) 1 - w_p(d)*Ts; % Euler
% a = @(d) exp(-w_p(d)*Ts); % Zero Order Hold
sigma_t = @(d) W.ti*d(1)*sqrt((1-a(d)^2)/(1-a(d))^2);
sigma_m = sqrt(Ts*W.q);
Q = @(d) diag([zeros(Lk-1,1); sigma_t(d)^2*w_p(d)^2]); % Covariance matrix of the process noise

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
y = zeros(Yk, N); % Initialize size of output vector for all k

% Generate the true state values
for k = 2:N
    xt(:,k) = f(xt(:,k-1),u(:,k-1),d(:,k-1)) + Ts*n(d(:,k-1));
    y(:,k-1) = h(xt(:,k-1),d(:,k-1)) + v(:,k-1);
end
for k=1:N
    p(k) = ve(xt(:,k),d(:,k))-xt(3,k);
    fr(k) = Fr(xt(:,k),d(:,k))/(B.B*B.m);
    tr(k) = (1-D.mu)*Tr(xt(:,k),d(:,k));
    %     pi(k) = vri(xt(:,k),1);
end
figure
plot(xt(1,:));
title("wr")
% figure
% plot(xt(2,:));
% title("xt")
% figure
% plot(xt(3,:));
% title("xtdot")
% figure
% plot(xt(6,:));
% title("xb")
figure
plot(xt(4,1:500));
title("yt")
figure
plot(xt(5,1:500));
title("ytdot")
figure
plot(xt(8,1:500));
title("yb")

figure
plot(1:N,y(2,:),1:N,y_me(2,:));
title("xb")
figure
plot(p(:))
title("vr")
figure
plot(1:N,3*xt(12,:)/(2*T.H*T.m),1:N,fr)
title("Tg & Fr")
% figure
% plot(1:N,3*xt(12,:)/(2*T.H*T.m),1:N,-(B.B*B.ky + T.k)*xt(4,:)/T.m - (B.B*B.cy + T.c)*xt(5,:)/T.m + B.B*B.ky*xt(8,:)/T.m + B.B*B.cy*xt(9,:)/T.m)
% title("Tg Applied")
% xlim([1 500])

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


function res = cp_ct(la,be,cl,lambdaVec,pitchVec)
[~,i_la] = min(abs(lambdaVec-abs(la)));
[~,i_be] = min(abs(pitchVec-be));
res = cl(i_la,i_be);
end