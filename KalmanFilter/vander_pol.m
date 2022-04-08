clc
clear all
close all

%% Vector Sizes

Ts = 0.05;
N = 150/0.05;

x_i = [2; 0]';

Lk = size(x_i,2); % Size of state vector
Yk = 1; % Size of measured vector
t = Ts*(1:N);

mu = 1;
f1 = @(x) x(2);
f2 = @(x) -x(1) - mu*x(2)*x(1)^2 + mu*x(2);

f = @(x) [f1(x); f2(x)]; % Nonlinear prediction
h = @(x) x(1);

Q = diag([0.02 0.1]);

R = 0.2;

n = sqrt(Q)*randn(Lk, N); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
% Initialize matrices
xt = zeros(Lk, N); % Initialize size of true state for all k
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k

% Runge-Kutta 4th order method
for k = 1:N-1
    k_1 = f(xt(:,k));
    k_2 = f(xt(:,k)+0.5*Ts*k_1);
    k_3 = f(xt(:,k)+0.5*Ts*k_2);
    k_4 = f(xt(:,k)+Ts*k_3);
    xt(:,k+1) = xt(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n(:,k);  % main equation
    
    yt(:,k) = h(xt(:,k)) + v(:,k);
end
yt(:,N) = h(xt(:,N)) + v(:,N);

%% Execute Unscented Kalman Filter
% Initialize state and covariance
xk = zeros(Lk, N); % Initialize size of state estimate for all k
xk(:,1) = x_i;
P0 = [0.01; 0.01].^2;
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
e = zeros(Yk,N);

%% Unscented Kalman Filter
P = P0; % Set first value of P to the initial P0
for k = 1:N-1
    % Step 1: Generate the sigma-points
    try
        sP = chol(P,'lower'); % Calculate square root of error covariance
    catch
        k
        break
    end
    
    % chi_p = "chi previous" = chi(k-1) % Untransformed sigma points
    chi_p = [xk(:,k), xk(:,k)*ones(1,Lk)+sqrt(Lk+lambda)*sP, ...
        xk(:,k)*ones(1,Lk)-sqrt(Lk+lambda)*sP]; % Untransformed sigma points
    
    % Step 2: Prediction Transformation
    % Propagate each sigma-point through prediction
    % chi_m = "chi minus" = chi(k|k-1)
    chi_m = zeros(Lk,n_sigma_p); % Transformed sigma points
    for j=1:n_sigma_p
        % Runge-Kutta 4th order method
        k_1 = f(chi_p(:,j));
        k_2 = f(chi_p(:,j)+0.5*Ts*k_1);
        k_3 = f(chi_p(:,j)+0.5*Ts*k_2);
        k_4 = f(chi_p(:,j)+Ts*k_3);
%         chi_m(:,j) = chi_p(:,j) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n(:,k);
        chi_m(:,j) = chi_p(:,j) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts;
%         chi_m(:,j) = chi_p(:,j) + Ts*f(chi_p(:,j),u_b(:,k));
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
%         psi_m(:,j) = h(chi_m(:,j)) + v(:,k);
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
    e(:,k+1) = yt(:,k+1) - y_m;
    xk(:,k+1) = x_m + K*e(:,k+1); % Update state estimate
    P = P_m - K*Pyy*K'; % Update covariance estimate
end

% xk = xk';
x_vl = {'$x_1$', '$x_2$'};
for i = 1:Lk
    figure
    subplot(1,2,1);
    %     plot(t,xt(i,:),'r-','LineWidth', 2);
    plot(t,xk(i,:),'b-', t,xt(i,:),'r-');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
%     ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF','True');
    title(['$Estimations $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);

    subplot(1,2,2);
    plot(t,xk(i,:)-xt(i,:),'b-');
    xlabel('$Time (s)$', 'Interpreter', 'latex', 'FontSize', 14);
%     ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF');
    title(['$Error $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end