function f = odeWindTurbine1(t,x,u,timeVector,To,Ac)
% StateFcnContinuous Evaluate the system ODEs

% To obtain the value of the inputs at the specified time
u = interp1(timeVector,u,t);

f(1) = -3*x(3)/(2*To.H*To.m) - To.c*x(1)/To.m - To.k*x(2)/To.m; % Tower edgewise acceleration
f(2) = x(1); % Tower edgewise velocity
f(3) = (u(1)-x(3))/Ac.tau; % Torque change in time

f = [f(1); f(2); f(3)];
end