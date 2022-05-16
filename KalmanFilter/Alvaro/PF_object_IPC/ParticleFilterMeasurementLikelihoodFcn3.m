function likelihood = ParticleFilterMeasurementLikelihoodFcn3(particles,measurement,B,D,To,M)
% ParticleFilterMeasurementLikelihoodFcn Measurement likelihood function
%
% likelihood = ParticleFilterMeasurementLikelihoodFcn(particles, measurement)
%
% Inputs:
%    particles - NumberOfStates-by-NumberOfParticles matrix that holds 
%                the particles
%
% Outputs:
%    likelihood - A vector with NumberOfParticles elements whose n-th
%                 element is the likelihood of the n-th particle
%
% See also extendedKalmanFilter, unscentedKalmanFilter

%   Copyright 2017 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

% Validate the sensor measurement
numberOfMeasurements = size(measurement,2); % Expected number of measurements
validateattributes(measurement, {'double'}, {'vector', 'numel', numberOfMeasurements}, ...
    'ParticleFilterMeasurementLikelihoodFcn1', 'measurement');

% Get all measurement hypotheses from particles
% predictedMeasurement = particles(1,:);
% R = [M.sigma_acc].^2;
% R = diag(R);
% v = sqrt(R)*randn(numberOfMeasurements, size(particles,2));
predictedMeasurement = zeros(numberOfMeasurements, size(particles,2)); % Initialize size of output vector for all k
for k = 1:size(particles,2)
    predictedMeasurement(:,k) = MeasurementFcn3(particles(:,k),B,D,To);
end

% Assume the ratio of the error between predicted and actual measurements
% follow a Gaussian distribution with zero mean, variance 0.2
mu = 0; % mean
% sigma = 0.15 * eye(numberOfMeasurements); % variance
% R = [M.sigma_enc; M.sigma_acc; M.sigma_acc; M.sigma_root; M.sigma_root;...
%     M.sigma_root; M.sigma_root; M.sigma_root; M.sigma_root; M.sigma_pit;...
%     M.sigma_pit; M.sigma_pit; M.sigma_pow; M.sigma_vane; M.sigma_azim].^2;
% sigma = diag(R);
R = [0.1; 0.01; 0.01; 0.01; 0.01;...
    0.01; 0.1; 0.1; 0.1; 0.01;...
    0.01; 0.01; 0.1; 0.1; 1].^2;
sigma = diag(R);

% Use multivariate Gaussian probability density function, calculate
% likelihood of each particle
numParticles = size(particles,2);
likelihood = zeros(numParticles,1);
C = det(2*pi*sigma) ^ (-0.5);
for kk=1:numParticles
    errorRatio = (predictedMeasurement(:,kk)-measurement')./predictedMeasurement(:,kk);
    v = errorRatio-mu;
    likelihood(kk) = C * exp(-0.5 * (v' / sigma * v) );
end
end