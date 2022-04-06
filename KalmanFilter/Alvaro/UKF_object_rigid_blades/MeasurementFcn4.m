function h = MeasurementFcn4(x,D,To,Ae,ct_l,lambdaVec,pitchVec)
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
ve = x(10) + x(9);
vr = ve - x(3);

lamb = (x(1)*Ae.Rr)/(vr);
ct = cp_ct(lamb,x(6),ct_l,lambdaVec,pitchVec);

Fr = 0.5*Ae.rho*Ae.Ar*(vr)^2*ct;

h = zeros(5,1);

h(1) = x(1);
h(2) = (Fr-To.c*x(3)-To.k*x(2))/To.m;
h(3) = (-(3/(2*To.H))*x(8)-To.c*x(5)-To.k*x(4))/To.m;
h(4) = D.eta*x(8)*x(1);
h(5) = vr;

h = [h(1); h(2); h(3); h(4); h(5)];

end