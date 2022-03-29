function f = StateFcnContinuous2Bladed(x,u,d,Ac,Ae,B,D,To,cp_l,ct_l,lambdaVec,pitchVec)
% StateFcnContinuous Evaluate the system ODEs

% w_p = v_m*pi/(2*W.L);
% ve = x(14) + x(13);
% vr = ve - x(3);

lamb = (x(1)*Ae.Rr-x(9))/(d(1)-x(7));
% lamb = @(x) (x(1)*Ae.Rr)/(vr(x));
cp = cp_ct(lamb,x(10),cp_l,lambdaVec,pitchVec);
ct = cp_ct(lamb,x(10),ct_l,lambdaVec,pitchVec);

Tr = x(8)*B.ky*2*B.l/3;
% Tr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x)-x(7))^3*cp(x)/x(1);
Fx = 0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^2*ct;
Fy = (0.5*Ae.rho*Ae.Ar*(d(1)-x(7))^3*cp*3)/(2*x(1)*B.l);

f = zeros(14,1);
%% Drive train
f(1) = (1-D.mu)*Tr/(D.Jr+D.Jg) - x(12)/(D.Jr+D.Jg);

%% Tower
f(2) = x(3); % Tower foreafter velocity
f(3) = -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.B*B.kx*x(6)/To.m + B.B*B.cx*x(7)/To.m; % Tower foreafter acceleration

f(4) = x(5); % Tower edgewise velocity
f(5) = -3*x(12)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.B*B.ky*x(8)/To.m + B.B*B.cy*x(9)/To.m; % Tower edgewise acceleration

%% Blades
f(6) = x(7); % Blade foreafter velocity
f(7) = Fx/(B.B*B.m) + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(7)/B.m; % Blade foreafter acceleration

f(8) = x(9); % Blade edgewise velocity
f(9) = Fy/(B.B*B.m) + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(8)/B.m - B.cy*x(9)/B.m; % Blade edgewise acceleration

%% Actuators
f(10) = x(11); % Pitch velocity
f(11) = Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(11) - Ac.omega^2*x(10); % Pitch acceleration
f(12) = (u(2)-x(12))/Ac.tau; % Torque change in time
% f(12) = (u(2)*x(1)^2-x(12))/Ac.tau; % Torque change in time (CONTROLLED)

%% Wind
% f(13) = -w_p*x(13); % Wind turbulence acceleration
% f(14) = 0; % Mean wind acceleration

f = [f(1); f(2); f(3); f(4); f(5); f(6); f(7); ...
    f(8); f(9); f(10); f(11); f(12)];
end