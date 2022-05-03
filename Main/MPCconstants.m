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
c4 = (B.B*c2 + To.k/To.m);
c5 = (B.B*c3 + To.c/To.m);

% xb1_dot
ws_ts = @(x,i) (To.r^2*(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2)^2 +...
    ((Ae.Rr*cos(x(27)+2*pi*i/3)+To.H)/To.H)^W.alpha); % Wind Share and Tower Shadow
vri_eq = @(x,i) x(26)*ws_ts(x,i) + x(25) - x(3)- x(9+i);
lambi_eq = @(x,i) (x(1)*Ae.Rr-x(15+i))/(vri_eq(x,i));
cpi_eq = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),cp_l,lambdaVec,pitchVec)/B.B;
cti_eq = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),ct_l,lambdaVec,pitchVec)/B.B;

% Gradients (dCt_dLambda, dCt_dPitch, dCp_dLambda, dCp_dPitch)
[dCp_dTheta_vec, dCp_dLambda_vec] = gradient(cp_l);
[dCt_dTheta_vec, dCt_dLambda_vec] = gradient(ct_l);
dCp_dTheta = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),dCp_dTheta_vec,lambdaVec,pitchVec)/B.B;
dCp_dLamb = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),dCp_dLambda_vec,lambdaVec,pitchVec)/B.B;
dCt_dTheta = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),dCt_dTheta_vec,lambdaVec,pitchVec)/B.B;
dCt_dLamb = @(x,i) cp_ct(lambi_eq(x,i),x(18+i),dCt_dLambda_vec,lambdaVec,pitchVec)/B.B;

dvri_dazim = @(x,i) x(26)*(To.r^2*((2*Ae.Rr^2*cos(x(27)+2*pi*i/3)*sin(x(27)+2*pi*i/3))...
    *(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 + To.xh^2)^2 - (Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 - To.xh^2)...
    *2*(2*Ae.Rr^2*cos(x(27)+2*pi*i/3)*sin(x(27)+2*pi*i/3))*(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 + To.xh^2))...
    /(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2 + To.xh^2)^4 ...
    + W.alpha*(-Ae.Rr*sin(x(27)+2*pi*i/3)/To.H)*((Ae.Rr*cos(x(27)+2*pi*i/3)+To.H)/To.H)^(W.alpha-1));

% lambda
r1 = @(x,i) Ae.Rr/vri_eq(x,i); % dlamb_dwr
r2 = @(x,i) -1/vri_eq(x,i); % dlamb_dydotbi
r3 = @(x,i) -(x(1)*Ae.Rr-x(15+i))*ws_ts(x,i)/vri_eq(x,i)^2; % dlamb_dvm
r4 = @(x,i) -(x(1)*Ae.Rr-x(15+i))*dvri_dazim(x,i)/vri_eq(x,i)^2; % dlamb_dazim
r5 = @(x,i) -(x(1)*Ae.Rr-x(15+i))/vri_eq(x,i)^2; % dlamb_dvt
r6 = @(x,i) (x(1)*Ae.Rr-x(15+i))/vri_eq(x,i)^2; % dlamb_dxdott
r7 = r6; % dlamb_dxdotbi


d1 = B.kx/B.m; % dx_bi
d2 = @(x,i) Ae.rho*Ae.Ar*(dCt_dLamb(x,i)*r7(x,i)*vri_eq(x,i)^2 - 2*cti_eq(x,i)*vri_eq(x,i))/(6*B.m) - B.cx/B.m; % dxdot_bi
d3 = d1; % dx_t
d4 = @(x,i) Ae.rho*Ae.Ar*(dCt_dLamb(x,i)*r6(x,i)*vri_eq(x,i)^2 - 2*cti_eq(x,i)*vri_eq(x,i))/(6*B.m) + B.cx/B.m; % dxdot_t
d5 = @(x,i) Ae.rho*Ae.Ar*(dCt_dLamb(x,i)*r3(x,i)*vri_eq(x,i)^2 + 2*cti_eq(x,i)*ws_ts(x,i)*vri_eq(x,i))/(6*B.m); % dv_m
d6 = @(x,i) Ae.rho*Ae.Ar*(dCt_dLamb(x,i)*r5(x,i)*vri_eq(x,i)^2 + 2*cti_eq(x,i)*vri_eq(x,i))/(6*B.m); % dv_t
d7 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^2*dCt_dTheta(x,i)/(6*B.m); % dpitch
d8 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^2*dCt_dLamb(x,i)*r2(x,i)/(6*B.m); % dydot_bi
d9 = @(x,i) Ae.rho*Ae.Ar*(2*cti_eq(x,i)*dvri_dazim(x,i)*vri_eq(x,i) + dCt_dLamb(x,i)*r4(x,i)*vri_eq(x,i)^2)/(6*B.m); % dazim
d10 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^2*dCt_dLamb(x,i)*r1(x,i)/(6*B.m); % dw_r

% yb1_dot
e1 = B.ky/B.m; % dy_bi
e2 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^3*(dCp_dLamb(x,i)*r2(x,i))/(4*B.m*x(1)*B.l) + B.cy/B.m; % dydot_bi
e3 = e1; % dy_t
e4 = B.cy/B.m; % dydot_t
e5 = @(x,i) Ae.rho*Ae.Ar*(dCp_dLamb(x,i)*r3(x,i)*vri_eq(x,i)^3 + 3*cpi_eq(x,i)*ws_ts(x,i)*vri_eq(x,i)^2)/(4*B.m*x(1)*B.l); % dv_m
e6 = @(x,i) Ae.rho*Ae.Ar*(dCp_dLamb(x,i)*r5(x,i)*vri_eq(x,i)^3 + 3*cpi_eq(x,i)*vri_eq(x,i)^2)/(4*B.m*x(1)*B.l); % dv_t
e7 = @(x,i) Ae.rho*Ae.Ar*(dCp_dLamb(x,i)*r6(x,i)*vri_eq(x,i)^3 - 3*cpi_eq(x,i)*vri_eq(x,i)^2)/(4*B.m*x(1)*B.l); % dxdot_t
e8 = @(x,i) Ae.rho*Ae.Ar*(dCp_dLamb(x,i)*r7(x,i)*vri_eq(x,i)^3 - 3*cpi_eq(x,i)*vri_eq(x,i)^2)/(4*B.m*x(1)*B.l); % dxdot_bi
e9 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^3*(dCp_dLamb(x,i)*r1(x,i)-cpi_eq(x,i)/x(1))/(4*B.m*x(1)*B.l); % dw_r
e10 = @(x,i) Ae.rho*Ae.Ar*(3*cpi_eq(x,i)*dvri_dazim(x,i)*vri_eq(x,i)^2 + dCp_dLamb(x,i)*r4(x,i)*vri_eq(x,i)^3)/(4*B.m*x(1)*B.l); % dazim
e11 = @(x,i) Ae.rho*Ae.Ar*vri_eq(x,i)^3*dCp_dTheta(x,i)/(4*B.m*x(1)*B.l); % dpitch

% Theta_dot
f1 = Ac.omega^2;
f2 = 2*Ac.omega*Ac.xi;

% Tg
g1 = 1/Ac.tau;
g2 = @(x,u) 2*u(4)*x(1)/Ac.tau;

% vt
h1 = @(x) W.w_p(x);

% My & Mx
p1 = B.kx*2*B.l/3;
p2 = B.ky*2*B.l/3;

% Pe
q1 = @(x) D.eta*x(1);
q2 = @(x) D.eta*x(24);


%%%%%%%%%%%%%%%
% CONSTRAINTS %
%%%%%%%%%%%%%%%

Sx = eye(Lk);
Sx(end-3,end-3) = S_means(1); %Tg

Sz = eye(Zk);
Sx(end,end) = S_means(2); %Pe

% Actuator Range Constraints
fc = [-1; 1];
f_rep = repmat({fc}, 1, Uk*Hu);
F_L = blkdiag(f_rep{:});
u_cons = [Ac.pitch_min; -Ac.pitch_max; Ac.pitch_min; -Ac.pitch_max;
    Ac.pitch_min; -Ac.pitch_max; Ac.Tg_min; -Ac.Tg_max];
U_L = u_cons;
for i=1:Hu-1
    U_L = [U_L; u_cons];
end
F = [F_L U_L];

% Actuator Slew Rate Constraints
ec = [-1; 1];
e_rep = repmat({ec}, 1, Uk-1);
e_rep{end+1} = zeros(2,1);
e_L = blkdiag(e_rep{:});

e_L = repmat({e_L}, 1, Hu);
E_L = blkdiag(e_L{:});
du_cons = [Ac.pitch_dot_min; -Ac.pitch_dot_max; Ac.pitch_dot_min; ...
    -Ac.pitch_dot_max; Ac.pitch_dot_min; -Ac.pitch_dot_max; zeros(2,1)];
dU_L = du_cons;
for i=1:Hu-1
    dU_L = [dU_L; du_cons];
end
E = [E_L dU_L];

% Constraints in Actuator Variables
g = [0; 0];
g_rep = repmat({g}, 1, Zk-7);
g_rep{1} = [-1;1];
for k=1:7
    g_rep{end+1} = [-1;1];
end
g_L = blkdiag(g_rep{:});

g_L = repmat({g_L}, 1, Hu);
G_L = blkdiag(g_L{:});
z_cons = cell2mat(struct2cell(Z_c));
Z_L = z_cons;
for i=1:Hu-1
    Z_L = [Z_L; z_cons];
end
G = [G_L Z_L];

ops = sdpsettings('solver','mosek','showprogress',0,'verbose',1,...
    'cachesolvers',1,'savedebug',0,'debug',1,'savesolverinput',0,'savesolveroutput',0);


