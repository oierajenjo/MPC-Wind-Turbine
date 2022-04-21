clc
clear all
close all

%% Obtain all variables
variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables
% run("variables.m")

% Generate inputs with controller
u_b = ones(2,N);
theta_f = 0;
[lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
u_b = [theta_f; K].*u_b;

%% Before filter execution

Lk = size(x_i_rigid,1); % Size of state vector
Yk = size(y_me_rigid,1); % Size of measured vector
Uk = size(u_b,1); % Size of imput vector

% Step 2: Define noise assumptions
w_p = @(x) x(10)*pi/(2*W.L);
ve = @(x) x(9) + x(10);
vr = @(x) ve(x) - x(3);

lamb = @(x) (x(1)*Ae.Rr)/(vr(x));
cp = @(x) cp_ct(lamb(x),x(6),cp_l,lambdaVec,pitchVec);
ct = @(x) cp_ct(lamb(x),x(6),ct_l,lambdaVec,pitchVec);

Tr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x))^3*cp(x)/x(1);
Fr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x))^2*ct(x);
%% Drive train
f1 = @(x) ((1-D.mu)*Tr(x) - x(8))/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3);
f3 = @(x) (Fr(x)-To.c*x(3)-To.k*x(2))/To.m;
f4 = @(x) x(5);
f5 = @(x) (-(3/(2*To.H))*x(8)-To.c*x(5)-To.k*x(4))/To.m;

%% Actuators
f6 = @(x) x(7); % Pitch velocity
f7 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(7) - Ac.omega^2*x(6); % Pitch acceleration
f8 = @(x,u) (u(2)*x(1)^2-x(8))/Ac.tau; % Torque change in time

%% Wind
f9 = @(x) -w_p(x)*x(9);
f10 = 0;

f = @(x,u) [f1(x); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x,u); f8(x,u); f9(x); f10]; % Nonlinear prediction

% h = @(x) (x);
h = @(x) [x(1); f3(x); f5(x); D.eta*x(8)*x(1); vr(x)];

rng(1);
a = @(x) 1- w_p(x)*Ts; % Euler
% a = @(x) exp(-(x(5)*pi/(2*L))*Ts); % Zero Order Hold
sigma_t = @(x) W.ti*x(10)*sqrt((1-a(x)^2)/(1-a(x))^2);
sigma_m = sqrt(W.q);
Q = @(x) diag([zeros(Lk-2,1); sigma_t(x)^2*w_p(x)^2; sigma_m^2]); % Covariance matrix of the process noise
temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_pow; 0.1].^2;
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
xt(:,1) = x_i_rigid; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k
% [xt,yt] = RK4(f,xt,h,yt,u_b,N,n,v,Ts);
for k = 1:N-1
    xt(:,k+1) = xt(:,k) + Ts*f(xt(:,k),u_b(:,k)) + Ts*n(xt(:,k));
    yt(:,k) = h(xt(:,k)) + v(:,k);
end
yt(:,N) = h(xt(:,N)) + v(:,N);

figure
plot(t,yt(1,:)',t,y_me_rigid(1,:));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('wr');

figure
plot(t,yt(2,:)',t,y_me_rigid(2,:));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('xtddot');

figure
plot(t,yt(3,:)',t,y_me_rigid(3,:));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('ytddot');

figure
plot(t,yt(4,:)',t,y_me_rigid(4,:));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('Pe');

figure
plot(t,yt(5,:)',t,y_me_rigid(5,:));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('vr');

figure
plot(t,xt(2,:)',t,data.Data(:,224));
legend('xTrue','Bladed')
% ylim([-2.6 2.6]);
ylabel('xt');

figure
plot(t,xt(4,:)',t,data.Data(:,225));
legend('xTrue','Bladed')
% ylim([-2.6 2.6]);
ylabel('yt');

figure
plot(t,xt(8,:)',t,data.Data(:,20));
legend('xTrue','Bladed')
% ylim([-2.6 2.6]);
ylabel('Tg');

%% Unscented Kalman Filter
% Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i_rigid;
P0 = [M.sigma_enc; M.sigma_tdef; M.sigma_tvel; M.sigma_tdef; M.sigma_tvel;...
    M.sigma_pit; M.sigma_pitvel; M.sigma_pow; 0.1; M.sigma_vane].^2;
P0 = diag(P0);
% P0 = 0.01*eye(Lk,Lk);

[xk,P,e] = UKF(f,h,Q,R,xk,yt,u_b,Lk,Yk,N,P0,Ts,v,n);

%% Display results
figure
plot(t,xt(1,:)',t,xk(1,:));
legend('True','UKF estimate')
% ylim([-2.6 2.6]);
ylabel('wr');

figure
plot(t,xt(2,:)',t,xk(2,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('xt');

figure
plot(t,xt(3,:)',t,xk(3,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('xtdot');

figure
plot(t,xt(4,:)',t,xk(4,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('yt');

figure
plot(t,xt(5,:)',t,xk(5,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('ytdot');

figure
plot(t,xt(6,:)',t,xk(6,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('theta');

figure
plot(t,xt(7,:)',t,xk(7,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('theta_dot');

figure
plot(t,xt(8,:)',t,xk(8,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('Tg');

figure
plot(t,xt(9,:)',t,xk(9,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('vt');

figure
plot(t,xt(10,:)',t,xk(10,:));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('vm');
