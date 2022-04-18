clc
clear all
close all

%% Obtain all variables
variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

Lk = size(x_i_rigid,1); % Size of state vector
Yk = size(y_me_rigid,1); % Size of measured vector
Uk = size(u_b,1); % Size of imput vector

%% Filter construction
pf = particleFilter(@ParticleFilterStateFcn4,@ParticleFilterMeasurementLikelihoodFcn4);
% Initialize it with 1000 particles around the mean x_i with 0.1 covariance.
n_part = 100;
initialize(pf, n_part, x_i_rigid, 0.01*eye(Lk));

% The filter uses as default option a state estimation method of 'mean',
% selected by the 'StateEstimationMethod' property and 'multinomial'
% resampling option via the 'ResamplingMethod' property, with certain
% characteristics seen in the 'ResamplingPolicy' property. Change if desired.

%% Start the estimation loop.
% This represents measurements arriving over time, step by step.

% Generate inputs with controller
u_b = ones(2,N);
theta_f = 0;
[lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
u_b = [theta_f; K].*u_b;

% Simulate the system for 600 seconds with the filter sample
% time 0.05 [s] to generate the true states of the system.
timeVector = 0:Ts:600-Ts;

rng(1); % Fix the random number generator for reproducible results
v_m = 6; % Constant mean wind speed
w_p = v_m*pi/(2*W.L); % Kaimal spectrum frequency
a = 1 - w_p*Ts; % Turbulent wind filter parameter using Euler
% a = @(x) exp(-w_p(x)*Ts); % Turbulent wind filter parameter using Zero Order Hold
sigma_t = W.ti*v_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind
sigma_m = sqrt(W.q); % Standard deviation mean wind
Q = diag([zeros(Lk-2,1); sigma_t^2*w_p^2; sigma_m^2]);
R = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_pow; 0.1].^2;
R = diag(R); % Covariance matrix of the measurement noise
n = sqrt(Q)*randn(Lk, N); % Generate random process noise (from assumed Q)
% v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

% [~,xTrue] = ode45(@(t,x) odeWindTurbine2(t,x,u_b,n,timeVector,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m), timeVector, x_i);
% xTrue = xTrue';

% Initialize matrices
xTrue = zeros(Lk, N); % Initialize size of true state for all k
xTrue(:,1) = x_i_rigid; % Set true initial state
% Runge-Kutta 4th order method
for k = 1:N-1
    k_1 = StateFcnContinuous4(xTrue(:,k),u_b(:,k),Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_2 = StateFcnContinuous4(xTrue(:,k)+0.5*Ts*k_1,u_b(:,k)+0.5*Ts,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_3 = StateFcnContinuous4(xTrue(:,k)+0.5*Ts*k_2,u_b(:,k)+0.5*Ts,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_4 = StateFcnContinuous4(xTrue(:,k)+Ts*k_3,u_b(:,k)+Ts,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    xTrue(:,k+1) = xTrue(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n(:,k);  % main equation
end

% Generate the measurements
% yt = xt(:,:);
yTrue = zeros(Yk, N); % Initialize size of output vector for all k
for k = 1:N
    yTrue(:,k) = MeasurementFcn4(xTrue(:,k),D,To,Ae,ct_l,lambdaVec,pitchVec);
end
% yMeas = yt .* (1+sqrt(R)*randn(size(yt))); % sqrt(R): Standard deviation of noise
% yMeas = yTrue + sqrt(R)*randn(size(yTrue));
% yMeas = yMeas';
yMeas = y_me_rigid';

figure
plot(timeVector,yTrue(1,:)',timeVector,yMeas(:,1));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('wr');

figure
plot(timeVector,yTrue(2,:)',timeVector,yMeas(:,2));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('xtddot');

figure
plot(timeVector,yTrue(3,:)',timeVector,yMeas(:,3));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('ytddot');

figure
plot(timeVector,yTrue(4,:)',timeVector,yMeas(:,4));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('Pe');

figure
plot(timeVector,yTrue(5,:)',timeVector,yMeas(:,5));
legend('yTrue','yMeas')
% ylim([-2.6 2.6]);
ylabel('vr');

figure
plot(timeVector,xTrue(2,:)',timeVector,data.Data(:,224));
legend('xTrue','Bladed')
% ylim([-2.6 2.6]);
ylabel('xt');

figure
plot(timeVector,xTrue(4,:)',timeVector,data.Data(:,225));
legend('xTrue','Bladed')
% ylim([-2.6 2.6]);
ylabel('yt');

figure
plot(timeVector,xTrue(8,:)',timeVector,data.Data(:,20));
legend('xTrue','Bladed')
% ylim([-2.6 2.6]);
ylabel('Tg');

% Estimate
xCorrectedPF = zeros(N,Lk);
for k=1:size(xTrue,2)
    % Use measurement y[k] to correct the particles for time k
    xCorrectedPF(k,:) = correct(pf,yMeas(k,:),D,To,Ae,M,ct_l,lambdaVec,pitchVec); % Filter updates and stores Particles[k|k], Weights[k|k]
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
    predict(pf,(u_b(:,k).*ones(Uk,n_part)),Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m); % Filter updates and stores Particles[k+1|k]
end

%% Plot the state estimates from particle filter:
figure
plot(timeVector,xTrue(1,:)',timeVector,xCorrectedPF(:,1));
legend('True','PF estimate')
% ylim([-2.6 2.6]);
ylabel('wr');

figure
plot(timeVector,xTrue(2,:)',timeVector,xCorrectedPF(:,2));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('xt');

figure
plot(timeVector,xTrue(3,:)',timeVector,xCorrectedPF(:,3));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('xtdot');

figure
plot(timeVector,xTrue(4,:)',timeVector,xCorrectedPF(:,4));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('yt');

figure
plot(timeVector,xTrue(5,:)',timeVector,xCorrectedPF(:,5));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('ytdot');

figure
plot(timeVector,xTrue(6,:)',timeVector,xCorrectedPF(:,6));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('theta');

figure
plot(timeVector,xTrue(7,:)',timeVector,xCorrectedPF(:,7));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('theta_dot');

figure
plot(timeVector,xTrue(8,:)',timeVector,xCorrectedPF(:,8));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('Tg');

figure
plot(timeVector,xTrue(9,:)',timeVector,xCorrectedPF(:,9));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('vt');

figure
plot(timeVector,xTrue(10,:)',timeVector,xCorrectedPF(:,10));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('vm');