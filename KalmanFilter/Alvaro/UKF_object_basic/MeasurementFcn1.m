function y = MeasurementFcn1(x,To)
% MeasurementFcn Example measurement function for discrete
% time nonlinear state estimators.
%
% y = MeasurementFcn(x)
%
% Inputs:
%    x - x[k], states at time k
% Outputs:
%    y - y[k], measurements at time k
%
% See also extendedKalmanFilter, unscentedKalmanFilter

%   Copyright 2016 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

% yk = xk(1)*(1+vk);

h = -3*x(3)/(2*To.H*To.m) - To.c*x(1)/To.m - To.k*x(2)/To.m;
y = h;

% y = x;
end