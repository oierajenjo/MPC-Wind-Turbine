function h = MeasurementFcn3(x,B,D,To)
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

ve = x(26) + x(25);
vr = ve - x(3);

h = zeros(15,1);

h(1) = x(1);
h(2) = -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.kx*sum(x(6:8))/To.m + B.cx*sum(x(9:11))/To.m;
h(3) = -3*x(24)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.ky*sum(x(12:14))/To.m + B.cy*sum(x(15:17))/To.m ;
h(4) = x(6)*B.kx*2*B.l/3;
h(5) = x(7)*B.kx*2*B.l/3;
h(6) = x(8)*B.kx*2*B.l/3;
h(7) = x(12)*B.ky*2*B.l/3;
h(8) = x(13)*B.ky*2*B.l/3;
h(9) = x(14)*B.ky*2*B.l/3;
h(10) = x(18);
h(11) = x(19);
h(12) = x(20);
h(13) = D.eta*x(24)*x(1);
h(14) = vr;
h(15) = x(27);

h = [h(1); h(2); h(3); h(4); h(5); h(6); h(7); h(8); h(9); h(10); h(11); h(12); h(13); h(14); h(15)];

% y = x;
end