function h = MeasurementFcn2(x,B,D,To)
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
ve = x(14) + x(13);
vr = ve - x(3);

h = zeros(7,1);

h(1) = x(1);
h(2) = -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.B*B.kx*x(6)/To.m + B.B*B.cx*x(7)/To.m;
h(3) = -3*x(12)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.B*B.ky*x(8)/To.m + B.B*B.cy*x(9)/To.m;
h(4) = -x(6)*B.kx*2*B.l/3;
h(5) = -x(8)*B.ky*2*B.l/3;
h(6) = D.eta*x(12)*x(1);
h(7) = vr;

h = [h(1); h(2); h(3); h(4); h(5); h(6); h(7)];

end