clc
clear all
close all

%% Obtain all variables
Copy_of_Basic_variables_CPC

%% Before filter execution
% Step 1: Define UT Scaling parameters and weight vectors
Lk = size(x_i,1); % Size of state vector
Yk = size(y_me,1); % Size of measured vector
Uk = size(u_b,1); % Size of imput vector

f1 = @(x) -3*x(3)/(2*To.H*To.m) - To.c*x(1)/To.m - To.k*x(2)/To.m; % Tower edgewise acceleration
f2 = @(x) x(1); % Tower edgewise velocity
f3 = @(x,u) (u(1)-x(3))/Ac.tau; % Torque change in time


f = @(x,u) [f1(x); f2(x); f3(x,u)]; % Nonlinear prediction
h = @(x) f1(x);

Q = diag(zeros(Lk,1)); % Covariance matrix of the process noise

temp = [M.sigma_acc].^2;
R = diag(temp); % Covariance matrix of measurement noise

% Step 3: Initialize state and covariance
x = zeros(Lk, N); % Initialize size of state estimate for all k
% x(:,1) = [0]; % Set initial state estimate
x(:,1) = x_i;

% Simulation Only: Calculate true state trajectory for comparison
% Also calculate measurement vector
% Var(QX) = QVar(X)Q' = sigma^4 -> Var(sqrt(Q)X) = sqrt(Q)Var(X)sqrt(Q)' = sigma^2
n = sqrt(Q)*randn(Lk, 1); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)
%w = zeros(1,N);
%v = zeros(1,N);
xt = zeros(Lk, N); % Initialize size of true state for all k
% xt(:,1) = zeros(Lk,1) + sqrt(P0)*randn(Lk,1); % Set true initial state
xt(:,1) = x_i; % Set true initial state
yt = zeros(Yk, N); % Initialize size of output vector for all k
yt(:,1) = y_me(:,1);
% Generate the true state values
% for k = 1:N
%     p(:,k) = f(xt(:,k),u_b(:,k));
%     xt(:,k+1) = xt(:,k) + Ts*(f(xt(:,k+1),u_b(:,k+1))+f(xt(:,k),u_b(:,k)))/2;
% end

% % Backward Euler's method 
%  for k = 1:N
%      wprime = xt(:,k) + Ts*f(xt(:,k),u_b(:,k));
%      xt(:,k+1) = xt(:,k) + Ts*f(xt(:,k),wprime);
%  end


% Runge-Kutta 4th order method
for k = 1:N-1  
    k_1 = f(xt(:,k),u_b(:,k));
    k_2 = f(xt(:,k)+0.5*Ts*k_1,u_b(:,k)+0.5*Ts);
    k_3 = f(xt(:,k)+0.5*Ts*k_2,u_b(:,k)+0.5*Ts);
    k_4 = f(xt(:,k)+Ts*k_3,u_b(:,k)+Ts);
    xt(:,k+1) = xt(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts;  % main equation
    
    yt(:,k) = h(xt(:,k)) + v(:,k);
end
yt(:,N) = h(xt(:,N)) + v(:,N);

t = Ts*(1:N);
figure
plot(t,xt(2,:), t,data.Data(:,225));
title("yt")
legend(["Us" "Bladed"])
% xlim([1 50])
figure
plot(t,xt(1,:), t,data.Data(1,231));
title("ytdot")
legend(["Us" "Bladed"])

%% Tower and blades edgewise
% A2 = [-To.c/To.m -To.k/To.m -3/(2*To.H*To.m); 1 0 0; 0 0 -1/Ac.tau];
% B2 = [0;0;1/Ac.tau];
% C2 = eye(3);
% D2 = zeros(3,1);
% 
% sys=ss(A2,B2,C2,D2);
% isstable(sys)
% ts = data.Data(:,25);
% figure
% [y_2,ts_2,x_2] = lsim(sys,u_b',ts,x_i);
% plot(ts_2,y_2)
% xlabel('Time (sec)')
% ylabel('System response')
% legend('yt','yt_{dot}','Tg')
% 
% sys=ss(eye(Lk)+Ts*A2,Ts*B2,C2,D2);
% isstable(sys)
