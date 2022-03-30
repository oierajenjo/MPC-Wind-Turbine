function f = StateFcnContinuous3(x,u,Ac,Ae,B,D,To,W,cp_l,ct_l,lambdaVec,pitchVec,v_m)
% StateFcnContinuous Evaluate the system ODEs

w_p = v_m*pi/(2*W.L);
ve = x(26) + x(25);
vr = ve - x(3);

for i=0:2
ve(i+1) = x(26)*(To.r^2*(Ae.Rr^2*(sin(wrapTo2Pi(x(27)+2*pi*i/3)))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(wrapTo2Pi(x(27)+2*pi*i/3)))^2)^2 +...
    ((Ae.Rr*cos(wrapTo2Pi(x(27)+2*pi*i/3))+To.H)/To.H)^W.alpha) + x(25);
vr(i+1) = ve(i+1) - x(3);

% lamb = @(x) (x(1)*Ae.Rr-mean(x(15:17)))/(vr(x)-mean(x(9:11)));
lamb(i+1) = (x(1)*Ae.Rr-x(15+i))/(vr(i+1)-x(9+i));

% cp = @(x) cp_ct(lamb(x),mean(x(18:20)),cp_l,lambdaVec,pitchVec);
cp(i+1) = cp_ct(lamb(i+1),x(18+i),cp_l,lambdaVec,pitchVec)/B.B;
ct(i+1) = cp_ct(lamb(i+1),x(18+i),ct_l,lambdaVec,pitchVec)/B.B;

Tr = (-x(12)-x(13)-x(14))*B.ky*2*B.l/3;
% Tr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x)-mean(x(9:11)))^3*cp(x)/x(1);
% Tr = @(x) (0.5*Ae.rho*Ae.Ar*((vri(x,0)-x(9))^3*cpi(x,0)+(vri(x,1)-x(10))^3*cpi(x,1)+(vri(x,2)-x(11))^3*cpi(x,2))/x(1));
Fx(i+1) = 0.5*Ae.rho*Ae.Ar*(vr(i+1)-x(9+i))^2*ct(i+1); % Thrust coefficient
Fy(i+1) = (0.5*Ae.rho*Ae.Ar*(vr(i+1)-x(9+i))^3*cp(i+1)*3)/(2*x(1)*B.l);
end

%% Drive train
f(1) = (1-D.mu)*Tr/(D.Jr+D.Jg) - x(24)/(D.Jr+D.Jg);

%% Tower
f(2) = x(3); % Tower foreafter velocity
f(3) = -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.kx*sum(x(6:8))/To.m + B.cx*sum(x(9:11))/To.m; % Tower foreafter acceleration

f(4) = x(5); % Tower edgewise velocity
f(5) = -3*x(24)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.ky*sum(x(12:14))/To.m + B.cy*sum(x(15:17))/To.m ; % Tower edgewise acceleration

%% Blades
f(6) = x(9); % Blade 1 foreafter velocity
f(7) = x(10); % Blade 2 foreafter velocity
f(8) = x(11); % Blade 3 foreafter velocity

f(9) = Fx(1)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(9)/B.m; % Blade 1 foreafter acceleration
f(10) = Fx(2)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(7)/B.m - B.cx*x(10)/B.m; % Blade 2 foreafter acceleration
f(11) = Fx(3)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(8)/B.m - B.cx*x(11)/B.m; % Blade 3 foreafter acceleration

f(12) = x(15); % Blade 1 edgewise velocity
f(13) = x(16); % Blade 2 edgewise velocity
f(14) = x(17); % Blade 3 edgewise velocity

f(15) = -Fy(1)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(12)/B.m - B.cy*x(15)/B.m; % Blade 1 edgewise acceleration
f(16) = -Fy(2)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(13)/B.m - B.cy*x(16)/B.m; % Blade 2 edgewise acceleration
f(17) = -Fy(3)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(14)/B.m - B.cy*x(17)/B.m; % Blade 3 edgewise acceleration

%% Actuators
f(18) = x(21); % Pitch 1 velocity
f(19) = x(22); % Pitch 2 velocity
f(20) = x(23); % Pitch 3 velocity
f(21) = Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(21) - Ac.omega^2*x(18); % Pitch 1 acceleration
f(22) = Ac.omega^2*u(2) - 2*Ac.omega*Ac.xi*x(22) - Ac.omega^2*x(19); % Pitch 2 acceleration
f(23) = Ac.omega^2*u(3) - 2*Ac.omega*Ac.xi*x(23) - Ac.omega^2*x(20); % Pitch 3 acceleration

% f(24) = (u(4)-x(24))/Ac.tau; % Torque change in time
f(24) = (u(4)*x(1)^2-x(24))/Ac.tau; % Torque change in time

%% Wind
f(25) = -w_p*x(25); % Wind turbulence acceleration
f(26) = 0; % Mean wind acceleration

%% Azimuth
f(27) = x(1); % Azimuth velocity

f = [f(1); f(2); f(3); f(4); f(5); f(6); f(7);...
    f(8); f(9); f(10); f(11); f(12); f(13); f(14); f(15);...
    f(16); f(17); f(18); f(19); f(20); f(21); f(22);...
    f(23); f(24); f(25); f(26); f(27)];
end