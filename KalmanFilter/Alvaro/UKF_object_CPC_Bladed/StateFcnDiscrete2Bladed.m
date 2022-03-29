function [x] = StateFcnDiscrete2Bladed(x,u,d,Ac,Ae,B,D,To,cp_l,ct_l,lambdaVec,pitchVec)
% Example state transition function for discrete-time nonlinear state
% estimators.
%
% xk1 = StateFcn(xk)
%
% Inputs:
%    xk - States x[k]
%
% Outputs:
%    xk1 - Propagated states x[k+1]
%
% See also extendedKalmanFilter, unscentedKalmanFilter

%   Copyright 2016 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

% % Euler integration of continuous-time dynamics x'=f(x) with sample time dt
% Ts = 0.05; % [s] Sample time
% x = x + StateFcnContinuous(x)*Ts;

% Runge-Kutta 4th order method with sample time Ts
Ts = 0.05; % [s] Sample time

k_1 = StateFcnContinuous2Bladed(x,u,d,Ac,Ae,B,D,To,cp_l,ct_l,lambdaVec,pitchVec);
k_2 = StateFcnContinuous2Bladed(x+0.5*Ts*k_1,u+0.5*Ts,d+0.5*Ts,Ac,Ae,B,D,To,cp_l,ct_l,lambdaVec,pitchVec);
k_3 = StateFcnContinuous2Bladed(x+0.5*Ts*k_2,u+0.5*Ts,d+0.5*Ts,Ac,Ae,B,D,To,cp_l,ct_l,lambdaVec,pitchVec);
k_4 = StateFcnContinuous2Bladed(x+Ts*k_3,u+Ts,d+Ts,Ac,Ae,B,D,To,cp_l,ct_l,lambdaVec,pitchVec);
x = x + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts;  % main equation
end