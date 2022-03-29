clc
clear all
close all

%% Obtain all variables
variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

%% Filter construction
pf = particleFilter(@ParticleFilterStateFcn2,@ParticleFilterMeasurementLikelihoodFcn2);
% Initialize it with 1000 particles around the mean x_i with 0.1 covariance.
initialize(pf, 1000, x_i, 0.1*eye(Lk));

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
sigma_m = sqrt(Ts*W.q); % Standard deviation mean wind
Q = diag([zeros(Lk-2,1); sigma_t^2*w_p^2; sigma_m^2]);
n = sqrt(Q)*randn(Lk, N); % Generate random process noise (from assumed Q)
% v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

[~,xTrue] = ode45(@(t,x) odeWindTurbine2(t,x,u_b,n,timeVector,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m), timeVector, x_i);
xTrue = xTrue';

% Generate the measurements
% yt = xt(:,:);
yTrue = zeros(Yk, N); % Initialize size of output vector for all k
for k = 1:N
    yTrue(:,k) = MeasurementFcn2(xTrue(:,k),B,D,To);
end
% yMeas = yt .* (1+sqrt(R)*randn(size(yt))); % sqrt(R): Standard deviation of noise
yMeas = y_me';

% Estimate
xCorrectedPF = zeros(N,Lk);
for k=1:size(xTrue,1)
    % Use measurement y[k] to correct the particles for time k
    xCorrectedPF(k,:) = correct(pf,yMeas(k,:),B,D,To); % Filter updates and stores Particles[k|k], Weights[k|k]
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
    predict(pf,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m); % Filter updates and stores Particles[k+1|k]
end

%% Plot the state estimates from particle filter:
figure
plot(timeVector,xTrue(1,:)',timeVector,xCorrectedPF(:,1));
legend('True','PF estimate')
% ylim([-2.6 2.6]);
ylabel('x_1');

figure
plot(timeVector,xTrue(2,:)',timeVector,xCorrectedPF(:,2));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('x_2');

figure
plot(timeVector,xTrue(3,:)',timeVector,xCorrectedPF(:,3));
legend('True','PF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('x_3');