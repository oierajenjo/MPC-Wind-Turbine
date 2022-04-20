% wr
a1 = (1-D.mu)*2*B.ky*B.l/(3*(D.Jr+D.Jg));
a2 = 1/(D.Jr+D.Jg);

% xt_dot
b1 = B.kx/To.m;
b2 = B.cx/To.m;
b3 = (B.B*b1 + To.k/To.m);
b4 = (B.B*b2 + To.c/To.m);

% yt_dot
c1 = 3/(2*To.H*To.m);
c2 = B.ky/To.m;
c3 = B.cy/To.m;
c4 = (B.B*c1 + To.k/To.m);
c5 = (B.B*c2 + To.c/To.m);

% xb1_dot
ws_ts = @(x,i) (To.r^2*(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2)^2 +...
    ((Ae.Rr*cos(x(27)+2*pi*i/3)+To.H)/To.H)^W.alpha); % Wind Share and Tower Shadow
vri_eq = @(x,i) x(26)*ws_ts(x,i) + x(25) - x(3)- x(9+i);
lambi_eq = @(x,i) (x(1)*Ae.Rr-x(15+i))/(vri_eq(x,i));
cpi_eq = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),cp_l,lambdaVec,pitchVec)/B.B;
cti_eq = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),ct_l,lambdaVec,pitchVec)/B.B;
dct_dlamb = 1; % Gradient
dct_dpitch = 1; % Gradient

dvri_dazim = @(x,i) x(26)*(To.r^2*((2*Ae.Rr^2*cos(x(27)+2*pi*i/3)*sin(x(27)+2*pi*i/3))...
    *(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 + To.xh^2)^2 - (Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 - To.xh^2)...
    *2*(2*Ae.Rr^2*cos(x(27)+2*pi*i/3)*sin(x(27)+2*pi*i/3))*(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 + To.xh^2))...
    /(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 + To.xh^2)^4 ...
    + W.alpha*(-Ae.Rr*sin(x(27)+2*pi*i/3)/To.H)*((Ae.Rr*sin(x(27)+2*pi*i/3)+To.H)/To.H)^(W.alpha-1));

dlamb_dvt = @(x,i) -(x(1)*Ae.Rr-x(15+i))/vri_eq(x,i)^2;
dlamb_dxdott = @(x,i) (x(1)*Ae.Rr-x(15+i))/vri_eq(x,i)^2;
dlamb_dxdotbi = dlamb_dxdott;
dlamb_dvm = @(x,i) -(x(1)*Ae.Rr-x(15+i))*ws_ts(x,i)/vri_eq(x,i)^2;
dlamb_dydotbi = @(x,i) -1/vri_eq(x,i);
dlamb_dazim = @(x,i) -(x(1)*Ae.Rr-x(15+i))*dvri_dazim(x,i)/vri_eq(x,i)^2;
dlamb_dwr = @(x,i) Ae.Rr/vri_eq(x,i);

d1 = B.kx/B.m; % dx_bi
d2 = @(x,i) Ae.rho*Ae.Ar*(dct_dlamb*dlamb_dxdotbi(x,i)*vri_eq(x,i)^2 - 2*cti_eq(x,i)*vri_eq(x,i) - B.cx/B.m)/(6*B.m); % dxdot_bi
d3 = d1; % dx_t
d4 = @(x,i) Ae.rho*Ae.Ar*(dct_dlamb*dlamb_dxdott(x,i)*vri_eq(x,i)^2 - 2*cti_eq(x,i)*vri_eq(x,i) + B.cx/B.m)/(6*B.m); % dxdot_t
d5 = @(x,i) Ae.rho*Ae.Ar*(dct_dlamb*dlamb_dvm(x,i)*vri_eq(x,i)^2 + 2*cti_eq(x,i)*ws_ts(x,i)*vri_eq(x,i))/(6*B.m); % dv_m
d6 = @(x,i) Ae.rho*Ae.Ar*(dct_dlamb*dlamb_dvt(x,i)*vri_eq(x,i)^2 + 2*cti_eq(x,i)*vri_eq(x,i))/(6*B.m); % dv_t
d7 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^2*dct_dpitch/(6*B.m); % dpitch
d8 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^2*dct_dlamb*dlamb_dydotbi(x,i)/(6*B.m); % dydot_bi
d9 = @(x,i) Ae.rho*Ae.Ar*(2*cti_eq(x,i)*dvri_dazim(x,i)*vri_eq(x,i) + dct_dlamb*dlamb_dazim(x,i)*vri_eq(x,i)^2)/(6*B.m); % dazim
d10 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^2*dct_dlamb*dlamb_dwr(x,i)/(6*B.m); % dw_r

% yb1_dot
dcp_dlamb = 1; % Gradient
dcp_dpitch = 1; % Gradient

e1 = B.ky/B.m; % dy_bi
e2 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^3*(dcp_dlamb*dlamb_dydotbi(x,i)/(4*x(1)*B.l) + B.cy/B.m); % dydot_bi
e3 = e1; % dy_t
e4 = B.cy/B.m; % dydot_t
e5 = @(x,i) Ae.rho*Ae.Ar*(dcp_dlamb*dlamb_dvm(x,i)*vri_eq(x,i)^3 + 3*cpi_eq(x,i)*ws_ts(x,i)*vri_eq(x,i)^2)/(4*x(1)*B.l); % dv_m
e6 = @(x,i) Ae.rho*Ae.Ar*(dcp_dlamb*dlamb_dvt(x,i)*vri_eq(x,i)^3 + 3*cpi_eq(x,i)*vri_eq(x,i)^2)/(4*x(1)*B.l); % dv_t
e7 = @(x,i) Ae.rho*Ae.Ar*(dcp_dlamb*dlamb_dxdott(x,i)*vri_eq(x,i)^3 - 3*cpi_eq(x,i)*vri_eq(x,i)^2)/(4*x(1)*B.l); % dxdot_t
e8 = @(x,i) Ae.rho*Ae.Ar*(dcp_dlamb*dlamb_dxdotbi(x,i)*vri_eq(x,i)^3 - 3*cpi_eq(x,i)*vri_eq(x,i)^2)/(4*x(1)*B.l); % dxdot_bi
e9 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^3*(dcp_dlamb*dlamb_dwr(x,i)-cpi_eq(x,i)/x(1))/(4*x(1)*B.l); % dw_r
e10 = @(x,i) Ae.rho*Ae.Ar*(3*cpi_eq(x,i)*dvri_dazim(x,i)*vri_eq(x,i)^2 + dcp_dlamb*dlamb_dazim(x,i)*vri_eq(x,i)^3)/(4*x(1)*B.l); % dazim
e11 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^3*dcp_dpitch/(4*x(1)*B.l); % dpitch

% Theta_dot
f1 = 2*Ac.omega*Ac.xi;
f2 = Ac.omega;

% Tg
g1 = 1/Ac.tau;

% vt
h1 = @(x) W.w_p(x);

% My & Mx
p1 = B.kx*2*B.l/3;
p2 = B.ky*2*B.l/3;

% Pe
q1 = @(x) D.eta*x(24);
q2 = @(x) D.eta*x(1);