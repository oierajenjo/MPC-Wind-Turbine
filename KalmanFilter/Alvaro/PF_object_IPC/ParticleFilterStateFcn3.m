function particles = ParticleFilterStateFcn3(particles,u,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m) 
% ParticleFilterStateFcn Example state transition function for particle
%                           filter
%
% Discrete-time approximation to system ODEs. 
% Sample time is 0.05s.
%
% predictedParticles = ParticleFilterStateFcn(particles)
%
% Inputs:
%    particles - Particles at current time. Matrix with dimensions
%                [NumberOfStates NumberOfParticles] matrix
%
% Outputs:
%    predictedParticles - Predicted particles for the next time step
%
% See also particleFilter

%   Copyright 2017 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

[numberOfStates, numberOfParticles] = size(particles);

% Time-propagate each particle
% Runge-Kutta 4th order method with sample time Ts
Ts = 0.05; % [s] Sample time

for kk=1:numberOfParticles
    k_1 = StateFcnContinuous3(particles(:,kk),u(:,kk),Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_2 = StateFcnContinuous3(particles(:,kk)+0.5*Ts*k_1,u(:,kk)+0.5*Ts,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_3 = StateFcnContinuous3(particles(:,kk)+0.5*Ts*k_2,u(:,kk)+0.5*Ts,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    k_4 = StateFcnContinuous3(particles(:,kk)+Ts*k_3,u(:,kk)+Ts,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
    particles(:,kk) = particles(:,kk) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts;
end

% Add process noise
v_m = 6; % Constant mean wind speed
w_p = v_m*pi/(2*W.L); % Kaimal spectrum frequency
a = 1 - w_p*Ts; % Turbulent wind filter parameter using Euler
% a = @(x) exp(-w_p(x)*Ts); % Turbulent wind filter parameter using Zero Order Hold
sigma_t = W.ti*v_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind
sigma_m = sqrt(W.q); % Standard deviation mean wind
Q = diag([zeros(numberOfStates-3,1); sigma_t^2*w_p^2; sigma_m^2; 0]); % Covariance matrix of the process noise
n = sqrt(Q)*randn(size(particles));
particles = particles + n;

% Add extra noise to every state variable
% extraNoise = 0.025*eye(numberOfStates);
% particles = particles + extraNoise * randn(size(particles));
end