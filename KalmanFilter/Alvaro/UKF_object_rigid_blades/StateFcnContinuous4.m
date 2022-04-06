function f = StateFcnContinuous4(x,u,Ac,Ae,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m)
% StateFcnContinuous Evaluate the system ODEs

w_p = v_m*pi/(2*W.L);
ve = x(10) + x(9);
vr = ve - x(3);

lamb = (x(1)*Ae.Rr)/(vr);
cp = cp_ct(lamb,x(6),cp_l,lambdaVec,pitchVec);
ct = cp_ct(lamb,x(6),ct_l,lambdaVec,pitchVec);

Tr = 0.5*Ae.rho*Ae.Ar*(vr)^3*cp/x(1);
Fr = 0.5*Ae.rho*Ae.Ar*(vr)^2*ct;

f = zeros(10,1);
%% Drive train
f(1) = ((1-D.mu)*Tr - x(8))/(D.Jr+D.Jg);

%% Tower
f(2) = x(3);
f(3) = (Fr-To.c*x(3)-To.k*x(2))/To.m;
f(4) = x(5);
f(5) = (-(3/(2*To.H))*x(8)-To.c*x(5)-To.k*x(4))/To.m;

%% Actuators
f(6) = x(7); % Pitch velocity
f(7) = Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(7) - Ac.omega^2*x(6); % Pitch acceleration
f(8) = (u(2)*x(1)^2-x(8))/Ac.tau; % Torque change in time

%% Wind
f(9) = -w_p*x(9);
f(10) = 0;

f = [f(1); f(2); f(3); f(4); f(5); f(6); f(7); ...
    f(8); f(9); f(10)];
end