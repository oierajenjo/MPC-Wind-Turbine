clc
clear all
close all
rng(1);

variables_IPC
load('BladedFiles\performancemap_data.mat')
Constant_variables
MPCconstants
syms x [1 Lk]
syms f [1 Lk]
syms cp ct

df1 = diff(f1,x12)*x12 + diff(f1,x13)*x13 + diff(f1,x14)*x14 + diff(f1,x24)*x24;
% int(f1(x12,x13,x14,x24))

lambda(x) = x(26)*(To.r^2*(Ae.Rr^2*(sin(x(27)))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(x(27)))^2)^2 +...
    ((Ae.Rr*cos(x(27))+To.H)/To.H)^W.alpha) + x(25) - x(3) -x(9);

dlambda = diff(lambda,x26)*x26 + diff(lambda,x27)*x27 + diff(lambda,x25)*x25 ...
    + diff(lambda,x3)*x3 + diff(lambda,x9)*x9;


f9(x)= 0.5*Ae.rho*Ae.Ar*cp*(x(26)*(To.r^2*(Ae.Rr^2*(sin(x(27)))^2-To.xh^2) ...
    /(To.xh^2+Ae.Rr^2*(sin(x(27)))^2)^2 + ((Ae.Rr*cos(x(27))+To.H)/To.H)^W.alpha) ...
    + x(25) - x(3) - x(9))^2*1/B.m  + B.kx*x(2)/B.m + B.cx*x(3)/B.m ...
    - B.kx*x(6)/B.m - B.cx*x(9)/B.m; % Blade 1 foreafter acceleration

df9 = diff(f9,x26)*x26 + diff(f9,x27)*x27 + diff(f9,x25)*x25 + diff(f9,x3)*x3 + ...
    diff(f9,x9)*x9 + diff(f9,x2)*x2 + diff(f9,x6)*x6;


f15(x)= -0.5*Ae.rho*Ae.Ar*(x(26)*(To.r^2*(Ae.Rr^2*(sin(x(27)))^2-To.xh^2) ...
    /(To.xh^2+Ae.Rr^2*(sin(x(27)))^2)^2 + ((Ae.Rr*cos(x(27))+To.H)/To.H)^W.alpha) ...
    + x(25) - x(3)-x(9))^3*3/(2*x(1)*B.l)/B.m + B.ky*x(4)/B.m + ...
    B.cy*x(5)/B.m - B.ky*x(12)/B.m - B.cy*x(15)/B.m;

df15 = diff(f15,x26)*x26 + diff(f15,x27)*x27 + diff(f15,x25)*x25 + diff(f15,x3)*x3 + ...
    diff(f15,x9)*x9 + diff(f15,x1)*x1 + diff(f15,x5)*x5 + diff(f15,x12)*x12 + ...
    diff(f15,x15)*x15  + diff(f15,x4)*x4;

ws_ts = @(x,i) (To.r^2*(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2)^2 +...
    ((Ae.Rr*cos(x(27)+2*pi*i/3)+To.H)/To.H)^W.alpha); % Wind Share and Tower Shadow
vei = @(x,i) x(26)*ws_ts(x,i) + x(25);
vri = @(x,i) vei(x,i) - x(3);

% lamb = @(x) (x(1)*Ae.Rr-mean(x(15:17)))/(vr(x)-mean(x(9:11)));
lambi = @(x,i) (x(1)*Ae.Rr-x(15+i))/(vri(x,i)-x(9+i));

% cp = @(x) cp_ct(lamb(x),mean(x(18:20)),cp_l,lambdaVec,pitchVec);
cpi = @(x,i) cp_ct(lambi(x,i),x(18+i),cp_l,lambdaVec,pitchVec)/B.B;
cti = @(x,i) cp_ct(lambi(x,i),x(18+i),ct_l,lambdaVec,pitchVec)/B.B;

% Tr = @(x) 0.5*Ae.rho*Ae.Ar*(vr(x)-mean(x(9:11)))^3*cp(x)/x(1);
% Tr = @(x) (0.5*Ae.rho*Ae.Ar*((vri(x,0)-x(9))^3*cpi(x,0)+(vri(x,1)-x(10))^3*cpi(x,1)+(vri(x,2)-x(11))^3*cpi(x,2))/x(1));
Fxi = @(x,i) 0.5*Ae.rho*Ae.Ar*(vri(x,i)-x(9+i))^2*cti(x,i); % Thrust coefficient
Fyi = @(x,i) 0.5*Ae.rho*Ae.Ar*(vri(x,i)-x(9+i))^3*cpi(x,i)*3/(2*x(1)*B.l);

%% Drive train
f1(x) = (1-D.mu)*(-(x(12)+x(13)+x(14))*B.ky*2*B.l/3)/(D.Jr+D.Jg) - x(24)/(D.Jr+D.Jg);
df1 = diff(f1,x12)*x12 + diff(f1,x13)*x13 + diff(f1,x14)*x14 + diff(f1,x24)*x24;

%% Tower
f2 = @(x) x(3); % Tower foreafter velocity
f3 = @(x) -(B.B*B.kx + To.k)*x(2)/To.m - (B.B*B.cx + To.c)*x(3)/To.m + B.kx*sum(x(6:8))/To.m + B.cx*sum(x(9:11))/To.m; % Tower foreafter acceleration

f4 = @(x) x(5); % Tower edgewise velocity
f5 = @(x) -3*x(24)/(2*To.H*To.m) - (B.B*B.ky + To.k)*x(4)/To.m - (B.B*B.cy + To.c)*x(5)/To.m + B.ky*sum(x(12:14))/To.m + B.cy*sum(x(15:17))/To.m ; % Tower edgewise acceleration

%% Blades
f6 = @(x) x(9); % Blade 1 foreafter velocity
f7 = @(x) x(10); % Blade 2 foreafter velocity
f8 = @(x) x(11); % Blade 3 foreafter velocity

f9 = @(x) Fxi(x,0)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(9)/B.m; % Blade 1 foreafter acceleration
f10 = @(x) Fxi(x,1)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(7)/B.m - B.cx*x(10)/B.m; % Blade 2 foreafter acceleration
f11 = @(x) Fxi(x,2)/B.m + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(8)/B.m - B.cx*x(11)/B.m; % Blade 3 foreafter acceleration

f12 = @(x) x(15); % Blade 1 edgewise velocity
f13 = @(x) x(16); % Blade 2 edgewise velocity
f14 = @(x) x(17); % Blade 3 edgewise velocity

f15 = @(x) -Fyi(x,0)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(12)/B.m - B.cy*x(15)/B.m; % Blade 1 edgewise acceleration
f16 = @(x) -Fyi(x,1)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(13)/B.m - B.cy*x(16)/B.m; % Blade 2 edgewise acceleration
f17 = @(x) -Fyi(x,2)/B.m + B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(14)/B.m - B.cy*x(17)/B.m; % Blade 3 edgewise acceleration

%% Actuators
f18 = @(x) x(21); % Pitch 1 velocity
f19 = @(x) x(22); % Pitch 2 velocity
f20 = @(x) x(23); % Pitch 3 velocity
f21 = @(x,u) Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(21) - Ac.omega^2*x(18); % Pitch 1 acceleration
f22 = @(x,u) Ac.omega^2*u(2) - 2*Ac.omega*Ac.xi*x(22) - Ac.omega^2*x(19); % Pitch 2 acceleration
f23 = @(x,u) Ac.omega^2*u(3) - 2*Ac.omega*Ac.xi*x(23) - Ac.omega^2*x(20); % Pitch 3 acceleration

% f24 = @(x,u) (u(4)-x(24))/Ac.tau; % Torque change in time
f24 = @(x,u) (u(4)*x(1)^2-x(24))/Ac.tau; % Torque change in time

%% Wind
f25 = @(x) -W.w_p(x)*x(25); % Wind turbulence acceleration
f26 = 0; % Mean wind acceleration

%% Azimuth
f27 = @(x) x(1); % Azimuth velocity
