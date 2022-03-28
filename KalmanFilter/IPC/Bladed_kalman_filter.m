clc
clear all
close all

%% Obtain all variables
Bladed_variables_IPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

% u_b = ones(4,N);
% theta_f = 0;
% [lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
% K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
% u_b = [theta_f; theta_f; theta_f; K].*u_b;

%% State transition functions
vei = @(x,d,i) d(1)*(To.r^2*(Ae.Rr^2*(sin((x(25)+2*pi*i/3)))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin((x(25)+2*pi*i/3)))^2)^2 +...
    ((Ae.Rr*cos((x(25)+2*pi*i/3))+To.H)/To.H)^W.alpha);
vri = @(x,d,i) vei(x,d,i) - x(3);

% lamb = @(x,d) (x(1)*Ae.Rr-mean(x(15:17)))/(d(1)-mean(x(9:11)));
lambi = @(x,d,i) (x(1)*Ae.Rr-x(15+i))/(vri(x,d,i)-x(9+i));

% cp = @(x,d) cp_ct(lamb(x,d),mean(x(18:20)),cp_l,lambdaVec,pitchVec);
cpi = @(x,d,i) cp_ct(lambi(x,d,i),x(18+i),cp_l,lambdaVec,pitchVec)/B.B;
cti = @(x,d,i) cp_ct(lambi(x,d,i),x(18+i),ct_l,lambdaVec,pitchVec)/B.B;


Tr = @(x,d) (x(12)+x(13)+x(14))*B.ky*2*B.l/3;
% Tr = @(x,d) 0.5*Ae.rho*Ae.Ar*(d(1)-x(3)-mean(x(9:11)))^3*cp(x,d)/x(1);
% Tr = @(x,d) (0.5*Ae.rho*Ae.Ar*((vri(x,d,0)-x(9))^3*cpi(x,d,0) +...
%      (vri(x,d,1)-x(10))^3*cpi(x,d,1) + (vri(x,d,2)-x(11))^3*cpi(x,d,2)))/x(1);
Fxi = @(x,d,i) 0.5*Ae.rho*Ae.Ar*(vri(x,d,i)-x(9+i))^2*cti(x,d,i); % Thrust coefficient
Fyi = @(x,d,i) (0.5*Ae.rho*Ae.Ar*(vri(x,d,i)-x(9+i))^3*cpi(x,d,i)*3)/(2*x(1)*B.l);

%% Drive train
f1 = @(x,d) (1-D.mu)*Tr(x,d)/(D.Jr+D.Jg) - x(24)/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3); % Tower foreafter velocity
f3 = @(x) -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.kx*sum(x(6:8))/To.m + B.cx*sum(x(9:11))/To.m; % Tower foreafter acceleration
% f3 = @(x) -To.Mr*9.81*To.xh*To.H/To.Jt - (B.B*B.kx + To.k)*(x(2)/To.H)*To.H/To.Jt - (B.B*B.cx + To.c)*(x(3)/To.H)*To.H/To.Jt + (B.kx*sum(x(6:8))/(2*B.l/3))*To.H/To.Jt + (B.cx*sum(x(9:11))/(2*B.l/3))*To.H/To.Jt; % Tower foreafter acceleration

f4 = @(x) x(5); % Tower edgewise velocity
f5 = @(x) -3*x(24)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.ky*sum(x(12:14))/To.m + B.cy*sum(x(15:17))/To.m ; % Tower edgewise acceleration

%% Blades
f6 = @(x) x(9); % Blade 1 foreafter velocity
f7 = @(x) x(10); % Blade 2 foreafter velocity
f8 = @(x) x(11); % Blade 3 foreafter velocity

f9 = @(x,d) Fxi(x,d,0)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(9)/B.m; % Blade 1 foreafter acceleration
f10 = @(x,d) Fxi(x,d,1)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(7)/B.m - B.cx*x(10)/B.m; % Blade 2 foreafter acceleration
f11 = @(x,d) Fxi(x,d,2)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(8)/B.m - B.cx*x(11)/B.m; % Blade 3 foreafter acceleration

f12 = @(x) x(15); % Blade 1 edgewise velocity
f13 = @(x) x(16); % Blade 2 edgewise velocity
f14 = @(x) x(17); % Blade 3 edgewise velocity

f15 = @(x,d) -Fyi(x,d,0)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(12)/B.m - B.cy*x(15)/B.m; % Blade 1 edgewise acceleration
f16 = @(x,d) -Fyi(x,d,1)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(13)/B.m - B.cy*x(16)/B.m; % Blade 2 edgewise acceleration
f17 = @(x,d) -Fyi(x,d,2)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(14)/B.m - B.cy*x(17)/B.m; % Blade 3 edgewise acceleration

%% Actuators
f18 = @(x) x(21); % Pitch 1 velocity
f19 = @(x) x(22); % Pitch 2 velocity
f20 = @(x) x(23); % Pitch 3 velocity
f21 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(21) - Ac.omega^2*x(18); % Pitch 1 acceleration
f22 = @(x,u) Ac.omega^2*u(2) - 2*Ac.omega*Ac.xi*x(22) - Ac.omega^2*x(19); % Pitch 2 acceleration
f23 = @(x,u) Ac.omega^2*u(3) - 2*Ac.omega*Ac.xi*x(23) - Ac.omega^2*x(20); % Pitch 3 acceleration

f24 = @(x,u) (u(4)-x(24))/Ac.tau; % Torque change in time
% f24 = @(x,u) (u(4)*x(1)^2-x(24))/Ac.tau; % Torque change in time

%% Azimuth
f25 = @(x) x(1); % Azimuth velocity

f = @(x,u,d) [f1(x,d); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x);...
    f8(x); f9(x,d); f10(x,d); f11(x,d); f12(x); f13(x); f14(x); f15(x,d);...
    f16(x,d); f17(x,d); f18(x); f19(x); f20(x); f21(x,u); f22(x,u);...
    f23(x,u); f24(x,u); f25(x)]; % Nonlinear prediction

h = @(x,d) [x(1); f3(x); f5(x); -x(6)*B.kx*2*B.l/3; -x(7)*B.kx*2*B.l/3;
    -x(8)*B.kx*2*B.l/3; -x(12)*B.ky*2*B.l/3; -x(13)*B.ky*2*B.l/3; ...
    -x(14)*B.ky*2*B.l/3; D.eta*x(24)*x(1); d(1); x(25)];

Q = diag(zeros(Lk,1)); % Covariance matrix of the process noise

temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_root; M.sigma_root;...
    M.sigma_root; M.sigma_root; M.sigma_root; M.sigma_root; M.sigma_pow;...
    M.sigma_vane; M.sigma_azim].^2;
R = diag(temp); % Covariance matrix of measurement noise


% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = sqrt(Q)*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

% Initialize matrices
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k
[xt,yt] = BRK4(f,xt,h,yt,u_b,d_b,N,n,v,Ts);

for k=1:N
    p1(k) = vri(xt(:,k), d_b(:,k),1);
    p2(k) = vri(xt(:,k), d_b(:,k),2);
    p3(k) = vri(xt(:,k), d_b(:,k),3);
    tr(k) = (1-D.mu)*Tr(xt(:,k), d_b(:,k));
end

figure
plot(t,xt(1,:),t,y_me(1,:));
title("wr")
legend(["Us" "Bladed"])
% 
figure
plot(t,xt(2,:), t, data.Data(:,224));
title("xt")
legend(["Us" "Bladed"])
% 
% % figure
% % plot(xt(3,1:500));
% % title("xtdot")
% 
figure
plot(t,xt(6,:),t,data.Data(:,85));
title("xb1")
legend(["Us" "Bladed"])
% 
figure
plot(t,xt(4,:),t,data.Data(:,225));
title("yt")
legend(["Us" "Bladed"])
% 
% % figure
% % plot(xt(5,1:500));
% % title("ytdot")
% 
figure
plot(t,xt(12,:),t,data.Data(:,86));
title("yb1")
legend(["Us" "Bladed"])
% 
% figure
% plot(t,tr(:),t,xt(24,:))
% title("Tr & Tg")
% legend("Tr", "Tg")
% 
% figure
% plot(t,p1(:),t,p2(:),t,p3(:))
% title("vr1, vr2, vr3")
% 
figure
plot(t,yt(4,:),t,y_me(4,:));
title("My")
legend(["Us" "Bladed"])

figure
plot(t,yt(7,:),t,y_me(7,:));
title("Mx")
legend(["Us" "Bladed"])

%% Execute Unscented Kalman Filter
% Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i;
P0 = [M.sigma_enc; M.sigma_tdef; M.sigma_tvel; M.sigma_tdef; M.sigma_tvel;...
    M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel; M.sigma_bvel;...
    M.sigma_bvel; M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel;... 
    M.sigma_bvel; M.sigma_bvel; M.sigma_pit; M.sigma_pit; M.sigma_pit;...
    M.sigma_pitvel; M.sigma_pitvel; M.sigma_pitvel; M.sigma_pow; M.sigma_azim].^2;
P0 = diag(P0);

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

%% Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0
for k = 1:N-1
    % Step 1: Generate the sigma-points
    sP = chol(P,'lower'); % Calculate square root of error covariance
    % chi_p = "chi previous" = chi(k-1) % Untransformed sigma points
    chi_p = [xk(:,k), xk(:,k)*ones(1,Lk)+sqrt(Lk+lambda)*sP, ...
        xk(:,k)*ones(1,Lk)-sqrt(Lk+lambda)*sP]; % Untransformed sigma points
    
    % Step 2: Prediction Transformation
    % Propagate each sigma-point through prediction
    % chi_m = "chi minus" = chi(k|k-1)
    chi_m = zeros(Lk,n_sigma_p); % Transformed sigma points
    for j=1:n_sigma_p
        chi_m(:,j) = chi_p(:,j) + Ts*f(chi_p(:,j),u_b(:,k),d_b(:,k));
    end
    
    x_m = chi_m*wm; % Calculate mean of predicted state
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
    for j=1:n_sigma_p
        psi_m(:,j) = h(chi_m(:,j),d_b(:,k));
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
    xk(:,k+1) = x_m + K*(y_me(:,k+1) - y_m); % Update state estimate
    P = P_m - K*Pyy*K'; % Update covariance estimate
end

%% Display results
result_display(t,Lk,xk,xt,x_ul,x_vl)
