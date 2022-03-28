clc
clear all
close all

%% Obtain all variables
Simple_Basic_variables_CPC
load('BladedFiles\performancemap_data.mat')
Constant_variables

%% States, measurements and inputs size
Lk = size(x_i,1); % Size of state vector
Yk = size(y_me,1); % Size of measured vector
Uk = size(u_b,1); % Size of imput vector

%% Unscented Kalman Filter Construction
% Your initial state guess at time k, utilizing measurements up to time k-1: xhat[k|k-1]
initialStateGuess = x_i; % xhat[k|k-1]
% Construct the filter
ukf = unscentedKalmanFilter(...
    @StateFcnDiscrete1,... % State transition function
    @MeasurementFcn1,... % Measurement function
    initialStateGuess,...
    'HasAdditiveMeasurementNoise',true);

% Add measurement noise
R = [M.sigma_acc].^2;
R = diag(R); % Covariance matrix of the measurement noise
ukf.MeasurementNoise = R;

% Add process noise
Q = diag(zeros(Lk,1)); % Covariance matrix of the process noise
ukf.ProcessNoise = Q;

%% Estimation Using predict and correct Commands

% Simulate the system for 600 seconds with the filter sample
% time 0.05 [s] to generate the true states of the system.
timeVector = 0:Ts:600-Ts;
rng(1); % Fix the random number generator for reproducible results
n = sqrt(Q)*randn(Lk, N); % Generate random process noise (from assumed Q)
v = sqrt(R)*randn(Yk, N); % Generate random measurement noise (from assumed R)

[~,xt] = ode45(@(t,x) odeWindTurbine1(t,x,u_b,timeVector,To,Ac), timeVector, x_i);
xt = xt';

% Generate the measurements
% yt = xt(:,:);
yt = zeros(Yk, N); % Initialize size of output vector for all k
for k = 1:N
    yt(:,k) = MeasurementFcn1(xt(:,k),To) + v(:,k);
end
% yMeas = yt .* (1+sqrt(R)*randn(size(yt))); % sqrt(R): Standard deviation of noise
yMeas = y_me';

% Preallocate space for data to analyze later
xCorrectedUKF = zeros(N,Lk); % Corrected state estimates
PCorrected = zeros(N,Lk,Lk); % Corrected state estimation error covariances
e = zeros(N,Yk); % Residuals (or innovations)

% Perform online estimation of the states x using the correct and predict
% commands. Provide generated data to the filter one time step at a time.
for k=1:N
    % Let k denote the current time.
    %
    % Residuals (or innovations): Measured output - Predicted output
    e(k,:) = yMeas(k,:) - MeasurementFcn1(ukf.State,To)'; % ukf.State is x[k|k-1] at this point
    % Incorporate the measurements at time k into the state estimates by
    % using the "correct" command. This updates the State and StateCovariance
    % properties of the filter to contain x[k|k] and P[k|k]. These values
    % are also produced as the output of the "correct" command.
    [xCorrectedUKF(k,:), PCorrected(k,:,:)] = correct(ukf,yMeas(k,:),To);
    % Predict the states at next time step, k+1. This updates the State and
    % StateCovariance properties of the filter to contain x[k+1|k] and
    % P[k+1|k]. These will be utilized by the filter at the next time step.
    predict(ukf,u_b(:,k)',To,Ac);
end

%% Unscented Kalman Filter Results and Validation
figure
plot(timeVector,xt(1,:)',timeVector,xCorrectedUKF(:,1));
legend('True','UKF estimate')
% ylim([-2.6 2.6]);
xlabel('Time [s]');
ylabel('x_1');

figure
plot(timeVector,xt(2,:)',timeVector,xCorrectedUKF(:,2));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('x_2');

figure
plot(timeVector,xt(3,:)',timeVector,xCorrectedUKF(:,3));
legend('True','UKF estimate')
% ylim([-3 1.5]);
xlabel('Time [s]');
ylabel('x_3');