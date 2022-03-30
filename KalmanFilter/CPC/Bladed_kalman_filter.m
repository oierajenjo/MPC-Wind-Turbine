clc
clear all
close all

%% Obtain all variables
Bladed_variables_CPC
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
lamb = @(x,d) (x(1)*Ae.Rr-x(9))/(d(1)-x(7));
cp = @(x,d) cp_ct(lamb(x,d),x(10),cp_l,lambdaVec,pitchVec);
ct = @(x,d) cp_ct(lamb(x,d),x(10),ct_l,lambdaVec,pitchVec);

Tr = @(x,d) x(8)*B.ky*2*B.l/3;
% Tr = @(x,d) 0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^3*cp(x,d)/x(1);
Fx = @(x,d) 0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^2*ct(x,d);
Fy = @(x,d) (0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^3*cp(x,d)*3)/(2*x(1)*B.l);

%% Drive train
f1 = @(x,d) (1-D.mu)*Tr(x,d)/(D.Jr+D.Jg) - x(12)/(D.Jr+D.Jg);

%% Tower
f2 = @(x) x(3); % Tower foreaft velocity
f3 = @(x) -(B.B*B.kx + To.k)*(x(2)+0.6)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.B*B.kx*x(6)/To.m + B.B*B.cx*x(7)/To.m; % Tower foreaft acceleration

f4 = @(x) x(5); % Tower sideways velocity
f5 = @(x) -3*x(12)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.B*B.ky*x(8)/To.m + B.B*B.cy*x(9)/To.m; % Tower sideways acceleration

%% Blades
f6 = @(x) x(7); % Blade flapwise velocity
f7 = @(x,d) Fx(x,d)/(B.B*B.m) + B.kx*(x(2)+0.6)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(7)/B.m; % Blade flapwise acceleration

f8 = @(x) x(9); % Blade edgewise velocity
f9 = @(x,d) B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(8)/B.m - B.cy*x(9)/B.m; % Blade edgewise acceleration

%% Actuators
f10 = @(x) x(11); % Pitch velocity
f11 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(11) - Ac.omega^2*x(10); % Pitch acceleration
f12 = @(x,u) (u(2)-x(12))/Ac.tau; % Torque change in time
% f12 = @(x,u) (u(2)*x(1)^2-x(12))/Ac.tau; % Torque change in time (CONTROLLED)

f = @(x,u,d) [f1(x,d); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x,d);...
    f8(x); f9(x,d); f10(x); f11(x,u); f12(x,u)]; % Nonlinear prediction

h = @(x,d) [x(1); f3(x); f5(x); -x(6)*B.kx*2*B.l/3; -x(8)*B.ky*2*B.l/3;...
    D.eta*x(12)*x(1); d(1)];

Q = diag(zeros(Lk,1)); % Covariance matrix of the process noise

temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_root; M.sigma_root;...
    M.sigma_pow; M.sigma_vane].^2;
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
    frx(k) = Fx(xt(:,k),d_b(:,k))/(B.B*B.m);
    fry(k) = Fy(xt(:,k),d_b(:,k))/(B.B*B.m);
    tr(k) = (1-D.mu)*Tr(xt(:,k),d_b(:,k))/(D.Jr+D.Jg);
end

t = Ts*(1:N);

figure
plot(t,xt(1,:), t,y_me(1,:));
title("wr")
legend(["Us" "Bladed"])
% xlim([1 50])

figure
plot(t,xt(2,:), t, data.Data(:,224));
title("xt")
legend(["Us" "Bladed"])
% xlim([1 50])

figure
plot(t,xt(6,:),t,mean([data.Data(:,85) data.Data(:,91) data.Data(:,97)], 2));
title("xb")
% xlim([1 50])

figure
plot(t,xt(7,:));
title("xbdot")
% xlim([1 50])

figure
plot(t,xt(4,:), t,data.Data(:,225));
title("yt")
legend(["Us" "Bladed"])
% xlim([1 50])

figure
plot(t,xt(8,:), t,mean([data.Data(:,86) data.Data(:,92) data.Data(:,98)], 2));
title("yb")
% xlim([1 50])

figure
plot(t,d_b(1,:))
title("vr")
% xlim([1 50])

figure
plot(t,yt(7,:),t,y_me(7,:))
title("vr")
legend(["Us" "Bladed"])
% xlim([1 50])

figure
plot(t,yt(4,:),t,y_me(4,:));
title("My")
legend(["Us","Bladed"])
% xlim([1 50])

figure
plot(t,yt(5,:),t,y_me(5,:));
title("Mx")
legend(["Us","Bladed"])
% xlim([1 50])

% figure
% plot(t,yt(4,:)/(B.B*B.m*B.l),t,y_me(4,:)/(B.B*B.m*B.l));
% title("xbddot")
% legend(["Us","Bladed"])
% % xlim([1 50])

%% Execute Unscented Kalman Filter
% Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i;
P0 = [M.sigma_enc; M.sigma_tdef; M.sigma_tvel; M.sigma_tdef; M.sigma_tvel;...
    M.sigma_bdef; M.sigma_bvel; M.sigma_bdef; M.sigma_bvel;... 
    M.sigma_pit; M.sigma_pitvel; M.sigma_pow].^2;
P0 = diag(P0);
% P0 = 0.01*eye(Lk,Lk); % Set initial error covariance

[xk,P,e] = BUKF(f,h,Q,R,xk,y_me,u_b,d_b,Lk,Yk,N,P0,Ts);

%% Display results
result_display(t,Lk,xk,xt,x_ul,x_vl)
