clc
clear all
close all

%% Obtain all variables
variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

%% For CONTROLLED
% u_b = ones(2,N);
% 
% theta_f = 0;
% [lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
% K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
% u_b = [theta_f; K].*u_b;

%% Before filter execution
% Step 2: Define noise assumptions
w_p = @(x) x(14)*pi/(2*W.L);
ve = @(x) x(14) + x(13);
vr = @(x) ve(x) - x(3);

lamb = @(x) (x(1)*Ae.Rr-x(9))/(vr(x)-x(7));
% lamb = @(x) (x(1)*Ae.Rr)/(vr(x));
cp = @(x) cp_ct(lamb(x),x(10),cp_l,lambdaVec,pitchVec);
ct = @(x) cp_ct(lamb(x),x(10),ct_l,lambdaVec,pitchVec);

Tr = @(x,d) x(8)*B.ky*2*B.l/3;
% Tr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x)-x(7))^3*cp(x)/x(1);
Fx = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x)-x(7))^2*ct(x);
Fy = @(x) (0.5*Ae.rho*Ae.Ar*(vr(x)-x(7))^3*cp(x)*3)/(2*x(1)*B.l);

%% Drive train
f1 = @(x) (1-D.mu)*Tr(x)/(D.Jr+D.Jg) - x(12)/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3); % Tower foreafter velocity
f3 = @(x) -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.B*B.kx*x(6)/To.m + B.B*B.cx*x(7)/To.m; % Tower foreafter acceleration

f4 = @(x) x(5); % Tower edgewise velocity
f5 = @(x) -3*x(12)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.B*B.ky*x(8)/To.m + B.B*B.cy*x(9)/To.m; % Tower edgewise acceleration

%% Blades
f6 = @(x) x(7); % Blade foreafter velocity
f7 = @(x) Fx(x)/(B.B*B.m) + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(7)/B.m; % Blade foreafter acceleration

f8 = @(x) x(9); % Blade edgewise velocity
f9 = @(x) Fy(x)/(B.B*B.m) + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(8)/B.m - B.cy*x(9)/B.m; % Blade edgewise acceleration

%% Actuators BIEN
f10 = @(x) x(11); % Pitch velocity
f11 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(11) - Ac.omega^2*x(10); % Pitch acceleration
f12 = @(x,u) (u(2)-x(12))/Ac.tau; % Torque change in time
% f12 = @(x,u) (u(2)*x(1)^2-x(12))/Ac.tau; % Torque change in time (CONTROLLED)

%% Wind
f13 = @(x) -w_p(x)*x(13); % Wind turbulence acceleration
f14 = 0; % Mean wind acceleration

f = @(x,u) [f1(x); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x);...
    f8(x); f9(x); f10(x); f11(x,u); f12(x,u); f13(x);  f14]; % Nonlinear prediction

h = @(x) [x(1); f3(x); f5(x); -x(6)*B.kx*2*B.l/3; -x(8)*B.ky*2*B.l/3;...
    D.eta*x(12)*x(1); vr(x)];

a = @(x) 1 - w_p(x)*Ts; % Euler
% a = @(x) exp(-w_p(x)*Ts); % Zero Order Hold
sigma_t = @(x) W.ti*x(14)*sqrt((1-a(x)^2)/(1-a(x))^2);
sigma_m = sqrt(Ts*W.q);
Q = @(x) diag([zeros(Lk-2,1); sigma_t(x)^2*w_p(x)^2; sigma_m^2]); % Covariance matrix of the process noise

temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_root; M.sigma_root;...
    M.sigma_pow; M.sigma_vane].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(x) sqrt(Q(x))*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

% Initialize matrices
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k
[xt,yt] = RK4(f,xt,h,yt,u_b,N,n,v,Ts);

for k=1:N
    p(k) = vr(xt(:,k));
    fx(k) = Fx(xt(:,k))/(B.B*B.m);
    fy(k) = Fy(xt(:,k))/(B.B*B.m);
    tr(k) = (1-D.mu)*Tr(xt(:,k))/(D.Jr+D.Jg);
end

t = Ts*(1:N);
figure
plot(t,xt(1,:), t,y_me(1,:));
title("wr")
legend(["Us" "Bladed"])
% xlim([1 50])
% figure
% plot(t,xt(2,:), t,-data.Data(:,224));
% title("xt")
% legend(["Us" "Bladed"])
% % xlim([1 50])
% % figure
% % plot(xt(3,:));
% % title("xtdot")
% figure
% plot(t,xt(6,:));
% title("xb")
% % xlim([1 50])
%
% figure
% plot(t,xt(4,:), t,data.Data(:,225));
% title("yt")
% legend(["Us" "Bladed"])
% % xlim([1 50])
% % figure
% % plot(xt(5,:));
% % title("ytdot")
% figure
% plot(t,xt(8,:));
% title("yb")
% % xlim([1 50])

figure
plot(t,p(:),t,y_me(7,:))
title("ve")
legend(["Us" "Bladed"])
% xlim([1 50])

figure
plot(t,yt(4,:),t,y_me(4,:));
title("Mx")
legend(["Us" "Bladed"])
% xlim([1 50])

figure
plot(t,yt(5,:),t,y_me(5,:));
title("My")
legend(["Us" "Bladed"])
% xlim([1 50])

%% Execute Unscented Kalman Filter
% Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i;
P0 = [M.sigma_enc; M.sigma_tdef; M.sigma_tvel; M.sigma_tdef; M.sigma_tvel;...
    M.sigma_bdef; M.sigma_bvel; M.sigma_bdef; M.sigma_bvel;... 
    M.sigma_pit; M.sigma_pitvel; M.sigma_pow; M.sigma_vane; M.sigma_vane].^2;
P0 = diag(P0);
% P0 = 0.01*eye(Lk,Lk); % Set initial error covariance

xk = UKF(f,h,Q,R,xk,y_me,u_b,Lk,Yk,N,P0,Ts);

%% Display results
result_display(t,Lk,xk,xt,x_ul,x_vl)

% figure
% plot(t,vr(xk),'b-',t,vr(xt,d),'r-');
% xlabel('Time [s]', 'FontSize', 14);
% ylabel('Velocity [m/s]', 'FontSize', 14);
% grid on;
% legend('UKF', 'True');
% title('Effective wind speed [v_r]', 'FontSize', 14);
% set(gcf, 'PaperOrientation','landscape');
% saveas(figure(6),'Figures/Kalman_ve.pdf');