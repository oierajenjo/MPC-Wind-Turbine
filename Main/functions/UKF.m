function [xk,P,e] = UKF(f,h,Q,R,x,y,u,kal,P,Ts,v,n,P0)
n_sigma_p = kal.n_sigma_p;
lambda = kal.lambda;
wm = kal.wm;
wc = kal.wc;
Lk = size(x,1); % Size of state vector
Yk = size(y,1); % Size of measured vector

%% Unscented Kalman Filter
% Step 1: Generate the sigma-points
try
    sP = chol(P,'lower');% Calculate square root of error covariance
catch ME
    disp('Matrix is not symmetric positive definite');
%     P = P0;
%     sP = chol(P,'lower');% Calculate square root of error covariance
    return
end

% chi_p = "chi previous" = chi(k-1) % Untransformed sigma points
chi_p = [x, x*ones(1,Lk)+sqrt(Lk+lambda)*sP, ...
    x*ones(1,Lk)-sqrt(Lk+lambda)*sP]; % Untransformed sigma points

% Step 2: Prediction Transformation
% Propagate each sigma-point through prediction
% chi_m = "chi minus" = chi(k|k-1)
chi_m = zeros(Lk,n_sigma_p); % Transformed sigma points
for j=1:n_sigma_p
    % Runge-Kutta 4th order method
%     [chi_m(:,j),~] = RK4(f,chi_p(:,j),u,h,n(chi_p(:,j)),v,Ts);
    [chi_m(:,j),~] = RK4(f,chi_p(:,j),u,h,zeros(Lk,1),zeros(Yk,1),Ts);
end

x_m = chi_m*wm; % Calculate mean of predicted state
% Calculate covariance of predicted state
P_m = Q; % A priori covariance estimate
% P_m = 0; % A priori covariance estimate
for j = 1:n_sigma_p
    P_m = P_m + wc(j)*(chi_m(:,j) - x_m)*(chi_m(:,j) - x_m)';
end

% Step 3: Observation Transformation
% Propagate each sigma-point through observation
% Initial velocity will be considered as 0, as we need it for
% obtaining the acceleration
psi_m = zeros(Yk,n_sigma_p);
for j=1:kal.n_sigma_p
    psi_m(:,j) = h(chi_m(:,j)) + v;
end
y_m = psi_m*wm; % Calculate mean of predicted output

% Calculate covariance of predicted output
% and cross-covariance between state and output
Pyy = R;
% Pyy = 0;
%Pxy = zeros(L,2);
Pxy = 0;
for j = 1:n_sigma_p
    Pyy = Pyy + wc(j)*(psi_m(:,j) - y_m)*(psi_m(:,j) - y_m)';
    Pxy = Pxy + wc(j)*(chi_m(:,j) - x_m)*(psi_m(:,j) - y_m)';
end

% Step 4: Measurement Update
K = Pxy/Pyy; % Calculate Kalman gain
e = y - y_m;
xk = x_m + K*e; % Update state estimate
P = P_m - K*Pyy*K'; % Update covariance estimate

end