clc
clear all
close all

%% Obtain all variables
variables_IPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

u_b = ones(4,N);
theta_f = 0;
[lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
u_b = [theta_f; theta_f; theta_f; K].*u_b;

%% Before filter execution
% Step 2: Define noise assumptions
w_p = @(x) 6*pi/(2*W.L);
ve = @(x) x(26) + x(25);
vr = @(x) ve(x) - x(3);

vei = @(x,i) x(26)*(To.r^2*(Ae.Rr^2*(sin(wrapTo2Pi(x(27)+2*pi*i/3)))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(wrapTo2Pi(x(27)+2*pi*i/3)))^2)^2 +...
    ((Ae.Rr*cos(wrapTo2Pi(x(27)+2*pi*i/3))+To.H)/To.H)^W.alpha) + x(25);
vri = @(x,i) vei(x,i) - x(3);

% lamb = @(x) (x(1)*Ae.Rr-mean(x(15:17)))/(vr(x)-mean(x(9:11)));
lambi = @(x,i) (x(1)*Ae.Rr-x(15+i))/(vri(x,i)-x(9+i));

% cp = @(x) cp_ct(lamb(x),mean(x(18:20)),cp_l,lambdaVec,pitchVec);
cpi = @(x,i) cp_ct(lambi(x,i),x(18+i),cp_l,lambdaVec,pitchVec)/B.B;
cti = @(x,i) cp_ct(lambi(x,i),x(18+i),ct_l,lambdaVec,pitchVec)/B.B;

Tr = @(x,d) -(x(12)+x(13)+x(14))*B.ky*2*B.l/3;
% Tr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x)-mean(x(9:11)))^3*cp(x)/x(1);
% Tr = @(x) (0.5*Ae.rho*Ae.Ar*((vri(x,0)-x(9))^3*cpi(x,0)+(vri(x,1)-x(10))^3*cpi(x,1)+(vri(x,2)-x(11))^3*cpi(x,2))/x(1));
Fxi = @(x,i) 0.5*Ae.rho*Ae.Ar*(vri(x,i)-x(9+i))^2*cti(x,i); % Thrust coefficient
Fyi = @(x,i) (0.5*Ae.rho*Ae.Ar*(vri(x,i)-x(9+i))^3*cpi(x,i)*3)/(2*x(1)*B.l);

%% Drive train
f1 = @(x) (1-D.mu)*Tr(x)/(D.Jr+D.Jg) - x(24)/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3); % Tower foreafter velocity
f3 = @(x) -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.kx*sum(x(6:8))/To.m + B.cx*sum(x(9:11))/To.m; % Tower foreafter acceleration

f4 = @(x) x(5); % Tower edgewise velocity
f5 = @(x) -3*x(24)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.ky*sum(x(12:14))/To.m + B.cy*sum(x(15:17))/To.m ; % Tower edgewise acceleration

%% Blades
f6 = @(x) x(9); % Blade 1 foreafter velocity
f7 = @(x) x(10); % Blade 2 foreafter velocity
f8 = @(x) x(11); % Blade 3 foreafter velocity

f9 = @(x) Fxi(x,0)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(9)/B.m; % Blade 1 foreafter acceleration
f10 = @(x) Fxi(x,1)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(7)/B.m - B.cx*x(10)/B.m; % Blade 2 foreafter acceleration
f11 = @(x) Fxi(x,2)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(8)/B.m - B.cx*x(11)/B.m; % Blade 3 foreafter acceleration

f12 = @(x) x(15); % Blade 1 edgewise velocity
f13 = @(x) x(16); % Blade 2 edgewise velocity
f14 = @(x) x(17); % Blade 3 edgewise velocity

f15 = @(x) -Fyi(x,0)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(12)/B.m - B.cy*x(15)/B.m; % Blade 1 edgewise acceleration
f16 = @(x) -Fyi(x,1)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(13)/B.m - B.cy*x(16)/B.m; % Blade 2 edgewise acceleration
f17 = @(x) -Fyi(x,2)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(14)/B.m - B.cy*x(17)/B.m; % Blade 3 edgewise acceleration

%% Actuators
f18 = @(x) x(21); % Pitch 1 velocity
f19 = @(x) x(22); % Pitch 2 velocity
f20 = @(x) x(23); % Pitch 3 velocity
f21 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(21) - Ac.omega^2*x(18); % Pitch 1 acceleration
f22 = @(x,u) Ac.omega^2*u(2) - 2*Ac.omega*Ac.xi*x(22) - Ac.omega^2*x(19); % Pitch 2 acceleration
f23 = @(x,u) Ac.omega^2*u(3) - 2*Ac.omega*Ac.xi*x(23) - Ac.omega^2*x(20); % Pitch 3 acceleration

% f24 = @(x,u) (u(4)-x(24))/Ac.tau; % Torque change in time
f24 = @(x,u) (u(4)*x(1)^2-x(24))/Ac.tau; % Torque change in time

%% Wind
f25 = @(x) -w_p(x)*x(25); % Wind turbulence acceleration
f26 = 0; % Mean wind acceleration

%% Azimuth
f27 = @(x) x(1); % Azimuth velocity

f = @(x,u) [f1(x); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x);...
    f8(x); f9(x); f10(x); f11(x); f12(x); f13(x); f14(x); f15(x);...
    f16(x); f17(x); f18(x); f19(x); f20(x); f21(x,u); f22(x,u);...
    f23(x,u); f24(x,u); f25(x); f26; f27(x)]; % Nonlinear prediction

h = @(x) [x(1); f3(x); f5(x); x(6)*B.kx*2*B.l/3; x(7)*B.kx*2*B.l/3;
    x(8)*B.kx*2*B.l/3; x(12)*B.ky*2*B.l/3; x(13)*B.ky*2*B.l/3; ...
    x(14)*B.ky*2*B.l/3; x(18); x(19); x(20); D.eta*x(24)*x(1); vr(x); x(27)];

rng(1);
a = @(x) 1 - w_p(x)*Ts; % Euler
% a = @(x) exp(-(x(5)*pi/(2*L))*Ts); % Zero Order Hold
sigma_t = @(x) W.ti*6*sqrt((1-a(x)^2)/(1-a(x))^2);
sigma_m = sqrt(W.q);
Q = @(x) diag([zeros(Lk-3,1); sigma_t(x)^2*w_p(x)^2; sigma_m^2; 0]); % Covariance matrix of the process noise

temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_root; M.sigma_root;...
    M.sigma_root; M.sigma_root; M.sigma_root; M.sigma_root; M.sigma_pit;...
    M.sigma_pit; M.sigma_pit; M.sigma_pow; M.sigma_vane; M.sigma_azim].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Step 3: Initialize state and covariance
% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(x) sqrt(Q(x))*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
% v = zeros(Yk, N);

%% Runge-Kutta 4th order method
% Initialize matrices
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k
[xt,yt] = RK4(f,xt,h,yt,u_b,N,n,v,Ts);
% Euler
% for k = 1:N-1
%     xt(:,k+1) = xt(:,k) + Ts*f(xt(:,k),u_b(:,k)) + Ts*n(xt(:,k)); 
%     yt(:,k) = h(xt(:,k)) + v(:,k);
% end
% yt(:,N) = h(xt(:,N)) + v(:,N);

figure
plot(t,yt(1,:)',t,y_me(1,:));
legend('$\omega_r$ (sim.)','$\omega_r$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Rotor speed');
ylabel('$\omega_r$ [rad/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,yt(2,:)',t,y_me(2,:));
legend('$\ddot{x}_{t}$ (sim.)','$\ddot{x}_{t}$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower fore-aft acceleration');
ylabel('$\ddot{x}_{t}$ [m/$s^2$]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,yt(3,:)',t,y_me(3,:));
legend('$\ddot{y}_{t}$ (sim.)','$\ddot{y}_{t}$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower sidewards acceleration');
ylabel('$\ddot{y}_{t}$ [m/$s^2$]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,yt(4,:)',t,y_me(4,:));
legend('$M_{y_1}$ (sim.)','$M_{y_1}$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Flapwise root bending moment (blade 1)');
ylabel('$M_{y_1}$ [Nm]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,yt(7,:)',t,y_me(7,:));
legend('$M_{x_1}$ (sim.)','$M_{x_1}$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Edgewise root bending moment (blade 1)');
ylabel('$M_{x_1}$ [Nm]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,yt(13,:)',t,y_me(13,:));
legend('$P_e$ (sim.)','$P_e$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Electrical power');
ylabel('$P_e$ [W]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,yt(14,:)',t,y_me(14,:));
legend('$v_r$ (sim.)','$v_r$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Relative wind speed at the rotor');
ylabel('$v_r$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(2,:)',t,data.Data(:,224));
legend('$x_t$ (sim.)','$x_t$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower fore-aft deflection');
ylabel('$x_t$ [m]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(4,:)',t,data.Data(:,225));
legend('$y_t$ (sim.)','$y_t$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Tower sidewards deflection');
ylabel('$y_t$ [m]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(Ts*(1:400),xt(6,1:400)',Ts*(1:400),data.Data(1:400,85));
legend('$x_{b_1}$ (sim.)','$x_{b_1}$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Flapwise deflection (blade 1)');
ylabel('$x_{b_1}$ [m]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

% figure
% plot(t,xt(6,:)',t,xt(7,:)',t,xt(8,:)');
% legend('$x_{b_1}$ (sim.)','$x_{b_2}$ (sim.)','$x_{b_3}$ (sim.)','Interpreter','latex')
% % ylim([-2.6 2.6]);
% title('Flapwise deflection (blade 1, 2 & 3)');
% ylabel('$x_{b}$ [m]','Interpreter','latex')
% xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(9,:)');
legend('$\dot{x}_{b_1}$ (sim.)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Flapwise velocity (blade 1)');
ylabel('$\dot{x}_{b_1}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

% figure
% plot(Ts*(1:400),xt(12,1:400)',Ts*(1:400),data.Data(1:400,86));
% legend('$y_{b_1}$ (sim.)','$y_{b_1}$ (Bladed)','Interpreter','latex')
% % ylim([-2.6 2.6]);
% title('Edgewise deflection (blade 1) - Euler');
% ylabel('$y_{b_1}$ [m]','Interpreter','latex')
% xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(12,:)',t,data.Data(:,86));
legend('$y_{b_1}$ (sim.)','$y_{b_1}$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Edgewise deflection (blade 1)');
ylabel('$y_{b}$ [m]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

% figure
% plot(t,xt(12,:)',t,xt(13,:)',t,xt(14,:)');
% legend('$y_{b_1}$ (sim.)','$y_{b_2}$ (sim.)','$y_{b_3}$ (sim.)','Interpreter','latex')
% % ylim([-2.6 2.6]);
% title('Edgewise deflection (blade 1, 2 & 3)');
% ylabel('$y_{b}$ [m]','Interpreter','latex')
% xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(15,:)');
legend('$\dot{y}_{b_1}$ (sim.)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Edgewise velocity (blade 1)');
ylabel('$\dot{y}_{b_1}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(18,:)',t,data.Data(:,34));
legend('$\theta_1 (sim.)$','$\theta_1$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Pitch angle (blade 1)');
ylabel('$\theta_1$ [rad]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(22,:)',t,data.Data(:,37));
legend('$\dot{\theta}_1$ (sim.)','$\dot{\theta}_1$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Pitch rate (blade 1)');
ylabel('$\dot{\theta}_1$ [rad/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

figure
plot(t,xt(24,:)',t,data.Data(:,20));
legend('$T_g$ (sim.)','$T_g$ (Bladed)','Interpreter','latex')
% ylim([-2.6 2.6]);
title('Generator torque');
ylabel('$T_g$ [Nm]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')

for k=1:N
    tr(k) = Tr(xt(:,k));
end
for k=1:N
    p1(k) = vri(xt(:,k),1);
    p2(k) = vri(xt(:,k),2);
    p3(k) = vri(xt(:,k),3);
end
figure
plot(Ts*(1:600),p1(1:600),Ts*(1:600),p2(1:600),Ts*(1:600),p3(1:600))
legend('$v_{r_1}$ (sim.)','$v_{r_2}$ (sim.)','$v_{r_3}$ (sim.)','Interpreter','latex')
title("Individual blade wind speed")
ylabel('$v_{r_i}$ [m/s]','Interpreter','latex')
xlabel('time [s]','Interpreter','latex')



figure
plot(t,tr/(D.Jr+D.Jg), t,xt(24,:)/(D.Jr+D.Jg))
legend('Tr','Tg')
xlim([0 10])

%% Unscented Kalman Filter
% Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i;
P0 = [M.sigma_enc; M.sigma_tdef; M.sigma_tvel; M.sigma_tdef; M.sigma_tvel;...
    M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel; M.sigma_bvel;...
    M.sigma_bvel; M.sigma_bdef; M.sigma_bdef; M.sigma_bdef; M.sigma_bvel;...
    M.sigma_bvel; M.sigma_bvel; M.sigma_pit; M.sigma_pit; M.sigma_pit;...
    M.sigma_pitvel; M.sigma_pitvel; M.sigma_pitvel; M.sigma_pow;...
    M.sigma_vane; 0.01; M.sigma_azim].^2;
P0 = diag(P0);
% P0 = 0.01*eye(Lk,Lk); 

[xk,P,e] = UKF(f,h,Q,R,xk,y_me,u_b,Lk,Yk,N,P0,Ts,v,n);

% % Construct the filter
% ukf = unscentedKalmanFilter(f,... % State transition function
%     h,... % Measurement function
%     x_i,...
%     'HasAdditiveMeasurementNoise',true);
% 
% ukf.MeasurementNoise = R;
% ukf.ProcessNoise = Q;
% % Preallocate space for data to analyze later
% xCorrectedUKF = zeros(N,Lk); % Corrected state estimates
% P = zeros(N,Lk,Lk); % Corrected state estimation error covariances
% P(1,:,:) = P0;
% e = zeros(Yk,N); % Residuals (or innovations)
% for k=1:N-1
%     % Let k denote the current time.
%     %
%     % Residuals (or innovations): Measured output - Predicted output
%     e(:,k) = y_me(:,k) - h(ukf.State); % ukf.State is x[k|k-1] at this point
%     % Incorporate the measurements at time k into the state estimates by
%     % using the "correct" command. This updates the State and StateCovariance
%     % properties of the filter to contain x[k|k] and P[k|k]. These values
%     % are also produced as the output of the "correct" command.
%     [xk(:,k+1), P(k+1,:,:)] = correct(ukf,y_me(:,k));
%     % Predict the states at next time step, k+1. This updates the State and
%     % StateCovariance properties of the filter to contain x[k+1|k] and
%     % P[k+1|k]. These will be utilized by the filter at the next time step.
%     predict(ukf,u_b(:,k));
% end


%% Display results
result_display(t,Lk,xk,xt,x_ul,x_vl)

% figure
% plot(t,vr(xk),'b-',t,vr(xt),'r-');
% xlabel('Time [s]', 'FontSize', 14);
% ylabel('Velocity [m/s]', 'FontSize', 14);
% grid on;
% legend('UKF', 'True');
% title('Effective wind speed [v_r]', 'FontSize', 14);
% set(gcf, 'PaperOrientation','landscape');
% saveas(figure(6),'Figures/Kalman_ve.pdf');