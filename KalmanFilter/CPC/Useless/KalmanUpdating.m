clc
clear all
close all

%% Obtain all variables
variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

%% Before filter execution
% Step 2: Define noise assumptions
lamb = @(x,d) (x(1)*Ae.Rr)/(d(1));
% lamb = @(x,d) (x(1)*Ae.Rr-x(9))/(d(1)-x(7));
cp = @(x,d) cp_ct(lamb(x,d),x(10),cp_l,lambdaVec,pitchVec);
% cpy = @(x,d) cp_ct(lamb(x,d),x(10),cp_l,lambdaVec,pitchVec);
ct = @(x,d) cp_ct(lamb(x,d),x(10),ct_l,lambdaVec,pitchVec);

Tr = @(x,d) 0.5*Ae.rho*Ae.Ar*(d(1))^3*cp(x,d)/x(1);
Fx = @(x,d) 0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^2*ct(x,d);
Fy = @(x,d) (0.5*Ae.rho*Ae.Ar*(d(1)-x(9))^3*cp(x,d)/x(1))/(2*Ae.Rr/3);

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
f12 = @(x,u) (u(2)-x(12))/Ac.tau; % Torque change in time


f = @(x,u,d) x + Ts*[f1(x,d); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x,d);...
    f8(x); f9(x,d); f10(x); f11(x,u); f12(x,u)]; % Nonlinear prediction

h = @(x,d) [x(1); f3(x); f5(x); B.l*B.m*f7(x,d); B.l*B.m*f9(x,d); ...
    D.eta*x(12)*x(1); d(1)];

Q = @(d) diag([zeros(Lk-1,1); 0]); % Covariance matrix of the process noise

temp = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_root; M.sigma_root;...
    M.sigma_pow; M.sigma_vane].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Step 3: Initialize state and covariance
% x = zeros(Lk, N); % Initialize size of state estimate for all k 
x = zeros(Lk, N); % Initialize size of state estimate for all k
x(:,1) = x_i;
% xt = zeros(Lk, N); % Initialize size of true state for all k
xt = x;
P0 = 0.01*eye(Lk,Lk); % Set initial error covariance

% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = @(d) sqrt(Q(d))*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

%% Execute Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0
for k = 2:N
    xt(:,k) = f(x(:,k-1),u_b(:,k-1),d_b(:,k-1)) + Ts*n(d_b(:,k-1));
    %yt(:,k) = h(xt(:,k),d_b(:,k)) + v(:,k);
    
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
    x(:,k) = x_m + K*(y_me(:,k) - y_m); % Update state estimate
    P = P_m - K*Pyy*K'; % Update covariance estimate
end

%% Display results
t = Ts*(1:k);
for i = 1:Lk
    figure(i)
    subplot(1,2,1);
    %     plot(t,xt(i,:),'r-','LineWidth', 2);
    plot(t,x(i,1:k),'b-', t,xt(i,1:k),'r-');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF','True');
    title(['$Estimations $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
    
    subplot(1,2,2);
    plot(t,x(i,1:k)-xt(i,1:k),'b-');
    xlabel('$Time (s)$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF');
    title(['$Error$' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end

% v_r = @(x,d) d(1,:)+x(13,:)-x(3,:);
% figure(6)
% plot(t,v_r(x,d_b),'b-',t,v_r(xt,d_b),'r-');
% xlabel('Time [s]', 'FontSize', 14);
% ylabel('Velocity [m/s]', 'FontSize', 14);
% grid on;
% legend('UKF', 'True');
% title('Effective wind speed [v_r]', 'FontSize', 14);
% set(gcf, 'PaperOrientation','landscape');
% saveas(figure(6),'Figures/Kalman_ve.pdf');

function res = cp_ct(la,be,cl,lambdaVec,pitchVec)
[~,i_la] = min(abs(lambdaVec-abs(la)));
[~,i_be] = min(abs(pitchVec-be));
res = cl(i_la,i_be);
end