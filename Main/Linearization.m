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

%% Drive train
f1(x) = (1-D.mu)*(-(x(12)+x(13)+x(14))*B.ky*2*B.l/3)/(D.Jr+D.Jg) - x(24)/(D.Jr+D.Jg);

df1 = diff(f1,x12)*x12 + diff(f1,x13)*x13 + diff(f1,x14)*x14 + diff(f1,x24)*x24;
% int(f1(x12,x13,x14,x24))

lambda(x) = x(26)*(To.r^2*(Ae.Rr^2*(sin(x(27)))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(x(27)))^2)^2 +...
    ((Ae.Rr*cos(x(27))+To.H)/To.H)^W.alpha) + x(25) - x(3) -x(9);

dlambda = diff(lambda,x26)*x26 + diff(lambda,x27)*x27 + diff(lambda,x25)*x25 ...
    + diff(lambda,x3)*x3 + diff(lambda,x9)*x9;


f9(x)= 0.5*Ae.rho*Ae.Ar*(x(26)*(To.r^2*(Ae.Rr^2*(sin(x(27)))^2-To.xh^2) ...
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

