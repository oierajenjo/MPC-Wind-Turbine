clc
clear all
close all

%% Obtain all variables
variables_IPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

%% Filter construction
pf = particleFilter(@ParticleFilterStateFcn3,@ParticleFilterMeasurementLikelihoodFcn3);
% Initialize it with 1000 particles around the mean x_i with 0.1 covariance.
n_part = 100;
initialize(pf, n_part, x_i, 0.1*eye(Lk));

% The filter uses as default option a state estimation method of 'mean',
% selected by the 'StateEstimationMethod' property and 'multinomial'
% resampling option via the 'ResamplingMethod' property, with certain
% characteristics seen in the 'ResamplingPolicy' property. Change if desired.

%% Start the estimation loop.
% This represents measurements arriving over time, step by step.

% Generate inputs with controller
u_b = ones(4,N);
theta_f = 0;
[lamb_opt, cp_opt] = cp_max(theta_f,cp_l,lambdaVec,pitchVec);
K = 0.5*Ae.rho*Ae.Rr^5*pi*cp_opt/lamb_opt^3;
u_b = [theta_f; theta_f; theta_f; K].*u_b;

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
Q = diag([zeros(Lk-3,1); sigma_t^2*w_p^2; sigma_m^2; 0]);
n = sqrt(Q)*randn(Lk, N); % Generate random process noise (from assumed Q)
% v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

% [~,xTrue] = ode45(@(t,x) odeWindTurbine2(t,x,u_b,n,timeVector,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m), timeVector, x_i);
% xTrue = xTrue';

% Initialize matrices
xTrue = zeros(Lk, N); % Initialize size of true state for all k
xTrue(:,1) = x_i; % Set true initial state
% Runge-Kutta 4th order method
for k = 1:N-1
    k_1 = StateFcnContinuous3(xTrue(:,k),u_b(:,k),Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_2 = StateFcnContinuous3(xTrue(:,k)+0.5*Ts*k_1,u_b(:,k)+0.5*Ts,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_3 = StateFcnContinuous3(xTrue(:,k)+0.5*Ts*k_2,u_b(:,k)+0.5*Ts,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_4 = StateFcnContinuous3(xTrue(:,k)+Ts*k_3,u_b(:,k)+Ts,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    xTrue(:,k+1) = xTrue(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n(:,k);  % main equation
end

% Generate the measurements
% yt = xt(:,:);
yTrue = zeros(Yk, N); % Initialize size of output vector for all k
for k = 1:N
    yTrue(:,k) = MeasurementFcn3(xTrue(:,k),B,D,To);
end
% yMeas = yt .* (1+sqrt(R)*randn(size(yt))); % sqrt(R): Standard deviation of noise
yMeas = y_me';

figure
plot(timeVector,yTrue(1,:)',timeVector,yMeas(:,1));
legend('yTrue','yMeas')
title('wr');
xlabel('s');
ylabel('wr');

figure
plot(timeVector,yTrue(2,:)',timeVector,yMeas(:,2));
legend('yTrue','yMeas')
title('xtddot');
xlabel('s');
ylabel('xtddot');

figure
plot(timeVector,yTrue(3,:)',timeVector,yMeas(:,3));
legend('yTrue','yMeas')
title('ytddot');
xlabel('s');
ylabel('ytddot');

figure
plot(timeVector,yTrue(4,:)',timeVector,yMeas(:,4));
legend('yTrue','yMeas')
title('My1');
xlabel('s');
ylabel('My1');

figure
plot(timeVector,yTrue(7,:)',timeVector,yMeas(:,7));
legend('yTrue','yMeas')
title('Mx1');
xlabel('s');
ylabel('Mx1');

figure
plot(timeVector,yTrue(10,:)',timeVector,yMeas(:,10));
legend('yTrue','yMeas')
title('Pe');
xlabel('s');
ylabel('Pe');

figure
plot(timeVector,yTrue(11,:)',timeVector,yMeas(:,11));
legend('yTrue','yMeas')
title('vr');
xlabel('s');
ylabel('vr');

figure
plot(timeVector,xTrue(2,:)',timeVector,data.Data(:,224));
legend('xTrue','Bladed')
title('xt');
xlabel('s');
ylabel('xt');

figure
plot(timeVector,xTrue(4,:)',timeVector,data.Data(:,225));
legend('xTrue','Bladed')
title('yt');
xlabel('s');
ylabel('yt');

figure
plot(timeVector,xTrue(6,:)',timeVector,data.Data(:,85));
legend('xTrue','Bladed')
title('xb1');
xlabel('s');
ylabel('xb1');

figure
plot(timeVector,xTrue(12,:)',timeVector,data.Data(:,86));
legend('xTrue','Bladed')
% ylim([-2.6 2.6]);
ylabel('yb1');

figure
plot(timeVector,xTrue(24,:)',timeVector,data.Data(:,20));
legend('xTrue','Bladed')
title('Tg');
xlabel('s');
ylabel('Tg');

% Estimate
xCorrectedPF = zeros(N,Lk);
for k=1:size(xTrue,2)
    % Use measurement y[k] to correct the particles for time k
    xCorrectedPF(k,:) = correct(pf,yMeas(k,:),B,D,To,M); % Filter updates and stores Particles[k|k], Weights[k|k]
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
    predict(pf,(u_b(:,k).*ones(Uk,n_part)),Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m); % Filter updates and stores Particles[k+1|k]
end

%% Plot the state estimates from particle filter:
figure
plot(timeVector,xTrue(1,:)',timeVector,xCorrectedPF(:,1));
legend('True','PF estimate')
title('wr');
xlabel('Time [s]');
ylabel('wr');

figure
plot(timeVector,xTrue(2,:)',timeVector,xCorrectedPF(:,2));
legend('True','PF estimate')
title('xt');
xlabel('Time [s]');
ylabel('xt');

figure
plot(timeVector,xTrue(3,:)',timeVector,xCorrectedPF(:,3));
legend('True','PF estimate')
title('xtdot');
xlabel('Time [s]');
ylabel('xtdot');

figure
plot(timeVector,xTrue(4,:)',timeVector,xCorrectedPF(:,4));
legend('True','PF estimate')
title('yt');
xlabel('Time [s]');
ylabel('yt');

figure
plot(timeVector,xTrue(5,:)',timeVector,xCorrectedPF(:,5));
legend('True','PF estimate')
title('ytdot');
xlabel('Time [s]');
ylabel('ytdot');

figure
plot(timeVector,xTrue(6,:)',timeVector,xCorrectedPF(:,6));
legend('True','PF estimate')
title('xb1');
xlabel('Time [s]');
ylabel('xb1');

figure
plot(timeVector,xTrue(9,:)',timeVector,xCorrectedPF(:,9));
legend('True','PF estimate')
title('xb1dot');
xlabel('Time [s]');
ylabel('xb1dot');

figure
plot(timeVector,xTrue(12,:)',timeVector,xCorrectedPF(:,12));
legend('True','PF estimate')
title('yb1');
xlabel('Time [s]');
ylabel('yb1');

figure
plot(timeVector,xTrue(15,:)',timeVector,xCorrectedPF(:,15));
legend('True','PF estimate')
title('yb1dot');
xlabel('Time [s]');
ylabel('yb1dot');

figure
plot(timeVector,xTrue(18,:)',timeVector,xCorrectedPF(:,18));
legend('True','PF estimate')
title('theta1');
xlabel('Time [s]');
ylabel('theta1');

figure
plot(timeVector,xTrue(21,:)',timeVector,xCorrectedPF(:,21));
legend('True','PF estimate')
title('theta_dot');
xlabel('Time [s]');
ylabel('theta1_dot');

figure
plot(timeVector,xTrue(24,:)',timeVector,xCorrectedPF(:,24));
legend('True','PF estimate')
title('Tg');
xlabel('Time [s]');
ylabel('Tg');

figure
plot(timeVector,xTrue(25,:)',timeVector,xCorrectedPF(:,25));
legend('True','PF estimate')
title('vt');
xlabel('Time [s]');
ylabel('vt');

figure
plot(timeVector,xTrue(26,:)',timeVector,xCorrectedPF(:,26));
legend('True','PF estimate')
title('vm');
xlabel('Time [s]');
ylabel('vm');