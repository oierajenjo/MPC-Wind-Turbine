function f = StateFcnContinuous1(x,u,To,Ac)
% StateFcnContinuous Evaluate the system ODEs

f(1) = -3*x(3)/(2*To.H*To.m) - To.c*x(1)/To.m - To.k*x(2)/To.m; % Tower edgewise acceleration
f(2) = x(1); % Tower edgewise velocity
f(3) = (u(1)-x(3))/Ac.tau; % Torque change in time

f = [f(1); f(2); f(3)];
end