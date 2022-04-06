function [x] = StateFcnDiscrete4(x,u,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m)
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

k_1 = StateFcnContinuous4(x,u,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
k_2 = StateFcnContinuous4(x+0.5*Ts*k_1,u+0.5*Ts,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
k_3 = StateFcnContinuous4(x+0.5*Ts*k_2,u+0.5*Ts,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
k_4 = StateFcnContinuous4(x+Ts*k_3,u+Ts,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m);
x = x + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts;  % main equation
end