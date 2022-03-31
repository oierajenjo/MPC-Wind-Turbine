clc
clear all
close all

%% Obtain all variables
Simple_Basic_variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

%% Filter construction
pf = particleFilter(@ParticleFilterStateFcn1,@ParticleFilterMeasurementLikelihoodFcn1);
% Initialize it with 1000 particles around the mean x_i with 0.1 covariance.
initialize(pf, 1000, x_i, 0.01*eye(Lk));

% The filter uses as default option a state estimation method of 'mean',
% selected by the 'StateEstimationMethod' property and 'multinomial'
% resampling option via the 'ResamplingMethod' property, with certain
% characteristics seen in the 'ResamplingPolicy' property. Change if desired.

%% Start the estimation loop.
% This represents measurements arriving over time, step by step.

% Simulate the system for 600 seconds with the filter sample
% time 0.05 [s] to generate the true states of the system.
timeVector = t;
% rng(1); % Fix the random number generator for reproducible results
% n = sqrt(Q)*randn(Lk, N); % Generate random process noise (from assumed Q)
% v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

[~,xTrue] = ode45(@(t,x) odeWindTurbine1(t,x,u_b,timeVector,To,Ac), timeVector, x_i);
xTrue = xTrue';

% Generate the measurements
% yTrue = xTrue(:,:);
yTrue = zeros(Yk, N); % Initialize size of output vector for all k
for k = 1:N
    yTrue(:,k) = MeasurementFcn1(xTrue(:,k),To);
end
% R = [0.01].^2;
% R = diag(R);
% yMeas = yTrue .* (1+sqrt(R)*randn(size(yTrue))); % sqrt(R): Standard deviation of noise
% yMeas = yMeas';
yMeas = y_me';

% Estimate
xCorrectedPF = zeros(N,Lk);
for k=1:size(xTrue,1)
    % Use measurement y[k] to correct the particles for time k
    xCorrectedPF(k,:) = correct(pf,yMeas(k,:),To); % Filter updates and stores Particles[k|k], Weights[k|k]
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
    predict(pf,(u_b(:,k).*ones(1000,Uk))',To,Ac); % Filter updates and stores Particles[k+1|k]
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

figure
plot(timeVector,yTrue(1,:)',timeVector,yMeas(:,1));
legend('yTrue','yMeas')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('y');