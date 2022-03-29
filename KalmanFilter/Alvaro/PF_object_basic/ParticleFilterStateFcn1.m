function particles = ParticleFilterStateFcn1(particles,u,To,Ac) 
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
    k_1 = StateFcnContinuous1(particles(:,kk),u(:,kk),To,Ac);
    k_2 = StateFcnContinuous1(particles(:,kk)+0.5*Ts*k_1,u(:,kk)+0.5*Ts,To,Ac);
    k_3 = StateFcnContinuous1(particles(:,kk)+0.5*Ts*k_2,u(:,kk)+0.5*Ts,To,Ac);
    k_4 = StateFcnContinuous1(particles(:,kk)+Ts*k_3,u(:,kk)+Ts,To,Ac);
    particles(:,kk) = particles(:,kk) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts;
end

% Add process noise
Q = diag(zeros(numberOfStates,1));
n = sqrt(Q)*randn(size(particles));
particles = particles + n;

% Add extra noise to every state variable
extraNoise = 0.025*eye(numberOfStates);
particles = particles + extraNoise * randn(size(particles));
end