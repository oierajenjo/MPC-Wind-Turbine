%% Drive train model constants
D.Jr = 321699000; % Rotor moment of inertia
D.Jg = 3.223e6; % Generator moment of inertia
% D.c = 0.005; % Drive train damping
% D.k = 1.409e10; % Drive train torsion stiffness
D.mu = 0.05; % Drive train mechanical losses (friction)
D.eta = 0.93; % Generator efficiency

var.D = D;

%% Blades model constants
B.r = 5.2/2;
B.l = 117.1836; % Blade length
B.Mb = 65566; % Blade mass
B.m = B.Mb/10+3*B.Mb*B.r^2/(20*B.l); % Blade equivalent mass 
B.d = 0.03; % Blade damping ratio
% B.fx = 0.541*sqrt(1-B.d^2); % Blade freq. flapwise
% B.fy = 0.636*sqrt(1-B.d^2); % Blade freq. edgewise
B.fx = 0.541; % Blade freq. flapwise
B.fy = 0.636; % Blade freq. edgewise
B.wx = 2*pi*B.fx; % Blade x direction natural pulsation
B.wy = 2*pi*B.fy; % Blade y direction natural pulsation
B.cx = B.d*2*B.m*B.wx; % Blade damping x direction
B.kx = B.wx^2*B.m; % Blade stiffness x direction
B.cy = B.d*2*B.m*B.wy; % Blade damping y direction
B.ky = B.wy^2*B.m; % Blade stiffness y direction
B.B = 3; % Blade amount
B.cc = 1/5; % Correction coefficient Fx
B.xd_min = -8.7;
B.xd_max = 10;
B.yd_min = -5.9;
B.yd_max = 6.1;

var.B = B;

%% Tower model constants
To.Mn = 630888; % Nacelle mass
To.Mt = 1086002; % Tower mass
To.Mr = 387198; % Rotor mass
To.m = To.Mn + To.Mt/3; % Tower equivalent mass
% To.m = 2475680-B.m*B.B; % Tower mass
To.d = 0.005; % Tower damping ratio
To.f = 0.18; % Tower freq. flapwise
To.c = To.d*2*To.m*2*pi*To.f; % Tower damping
To.k = (2*pi*To.f)^2*To.m; % Tower stiffness
To.Ht = 144.582; % Tower height
To.r_top = 3.25; % Tower top radius
To.r_base = 5; % Tower base radius
To.H = To.Ht + 4.34799; % Hub height
To.r = (To.r_top-To.r_base)*(To.H-B.l)/To.H + To.r_base; % Tower radius
To.xh = 10.93; % Hub overhang
To.Jt = To.Mn*To.H^2 + To.Mt*To.H^2/3;
rr = sqrt(To.xh^2+To.H^2);
rn = sqrt(3.945^2+To.H^2);
alpha_r = sin(acos(To.H/sqrt(To.xh^2+To.H^2)));
alpha_n = sin(acos(To.H/sqrt(3.945^2+To.H^2)));
To.xtoff = (rr*To.Mr*9.807*alpha_r+rn*To.Mn*9.807*alpha_n)/(To.k*(rr*0.37+rn*0.63)*sin(alpha_n*0.37+alpha_r*0.63));
To.xd_min = 0.20;
To.xd_max = 0.20;
To.yd_min = 0.15;
To.yd_max = 0.15;

var.To = To;

%% Aerodynamic model constants
Ae.rho = 1.225; % Density of the air
Ae.Rr = 241.996/2; % Rotor radius
Ae.Ar = pi*Ae.Rr^2; % Rotor area

var.Ae = Ae;

%% Wind model constants
Ts = 0.05; % Sampling time
W.ti = 0.15; % Turbulence intensity
W.q = 2^2/600; % Incremental variance mean wind speed
W.mu_m = 6; % Fixed mean wind speed: 10 m/s
W.L = 340.2;
W.alpha = 0.15; % Wind shear exponent for smooth terrain
W.mu_v = mean(data.Data(:,59));
W.w_p = @(x) W.mu_v*pi/(2*W.L);
W.a = @(x) 1 - W.w_p(x)*Ts; % Euler
W.sigma_t = @(x) W.ti*W.mu_v*sqrt((1-W.a(x)^2)/(1-W.a(x))^2);
W.sigma_m = sqrt(W.q);

var.W = W;

%% Actuator constants
Ac.omega = 2.4*pi; % Natural frequency of pitch actuator model
Ac.xi = 0.8; % Damping factor of pitch actuator model
Ac.tau = 0.1; % Generator time constant
Ac.pitch_min = 0;
Ac.pitch_max = pi/2;
Ac.Tg_min = 1;
Ac.Tg_max = 2.159e7;
Ac.pitch_dot = deg2rad(9); % Pitch angle max and min angular speed
Ac.omega_min = convangvel(5,'rpm', 'rad/s');
Ac.omega_opt = convangvel(7.56,'rpm', 'rad/s');
Ac.omega_max = convangvel(10,'rpm', 'rad/s');

Ac.Pe_min = D.eta*Ac.Tg_min*Ac.omega_min;
Ac.Pe_opt = 15*10^6;
Ac.Pe_max = D.eta*Ac.Tg_max*Ac.omega_max;

var.Ac = Ac;

%% Measurement constants
M.sigma_enc = 0.017;
M.sigma_acc = 0.04;
M.sigma_root = 0.01; % ¿?
M.sigma_pow = 0.035;
M.sigma_vane = 1;
M.sigma_azim = 0.01;

M.sigma_tdef = 0.01; % ¿?
M.sigma_tvel = 0.01; % ¿?
M.sigma_bdef = 0.01; % ¿?
M.sigma_bvel = 0.01; % ¿?
M.sigma_pit = 0.01; % ¿?
M.sigma_pitvel = 0.01; % ¿?

var.M = M;

%% Vector Sizes
Lk = size(x_i,1); % Size of state vector
Yk = size(y_me,1); % Size of measured vector
Uk = size(u_b,1); % Size of imput vector
t = Ts*(1:N);

%% Kalman variables
kal.alpha = 1; % Primary scaling parameter
kal.beta = 2; % Secondary scaling parameter (Gaussian assumption)
kal.kappa = 0; % Tertiary scaling parameter
kal.lambda = kal.alpha^2*(Lk+kal.kappa) - Lk;
kal.n_sigma_p = 2*Lk + 1; % Number of sigma points
kal.wm = ones(kal.n_sigma_p,1)*1/(2*(Lk+kal.lambda)); % Weight for transformed mean
kal.wc = kal.wm; % Weight for transformed covariance
kal.wm(1) = kal.lambda/(kal.lambda+Lk);
kal.wc(1) = kal.lambda/(kal.lambda+Lk) + 1 - kal.alpha^2 + kal.beta;

%% MPC variables
Hp = 10; % Prediction Horizon
Hu = 10; % Control Horizon
Hw = 1; % Window parameter

% Constraints limit values
Z_c.omega_min = Ac.omega_min;
Z_c.omega_max = -Ac.omega_max;
Z_c.xtd_min = To.xd_min;
Z_c.xtd_max = -To.xd_max;
Z_c.ytd_min = To.yd_min;
Z_c.ytd_max = -To.yd_max;

Z_c.xbid_min1 = B.xd_min;
Z_c.xbid_max1 = -B.xd_max;
Z_c.xbid_min2 = B.xd_min;
Z_c.xbid_max2 = -B.xd_max;
Z_c.xbid_min3 = B.xd_min;
Z_c.xbid_max3 = -B.xd_max;

Z_c.ybid_min1 = B.yd_min;
Z_c.ybid_max1 = -B.yd_max;
Z_c.ybid_min2 = B.yd_min;
Z_c.ybid_max2 = -B.yd_max;
Z_c.ybid_min3 = B.yd_min;
Z_c.ybid_max3 = -B.yd_max;

Z_c.pitchi_min1 = Ac.pitch_min;
Z_c.pitchi_max1 = -Ac.pitch_max;
Z_c.pitchi_min2 = Ac.pitch_min;
Z_c.pitchi_max2 = -Ac.pitch_max;
Z_c.pitchi_min3 = Ac.pitch_min;
Z_c.pitchi_max3 = -Ac.pitch_max;

Z_c.pitchid_min1 = Ac.pitch_dot;
Z_c.pitchid_max1 = -Ac.pitch_dot;
Z_c.pitchid_min2 = Ac.pitch_dot;
Z_c.pitchid_max2 = -Ac.pitch_dot;
Z_c.pitchid_min3 = Ac.pitch_dot;
Z_c.pitchid_max3 = -Ac.pitch_dot;

Z_c.Pe_min = Ac.Pe_min;
Z_c.Pe_max = -Ac.Pe_opt;

% Z_c.omega_r_min = convangvel(5,'rpm', 'rad/s');
% Z_c.omega_r_max = - convangvel(10,'rpm', 'rad/s');
% Z_c.xtdd_min = 0;
% Z_c.xtdd_max = -0;
% Z_c.ytdd_min = 0;
% Z_c.ytdd_max = -0;
% Z_c.My_min1 = B.m*B.l*B.xdd_min; Z_c.My_max1 = -B.m*B.l*B.xdd_max;
% Z_c.My_min2 = B.m*B.l*B.xdd_min; Z_c.My_max2 = -B.m*B.l*B.xdd_max;
% Z_c.My_min3 = B.m*B.l*B.xdd_min; Z_c.My_max3 = -B.m*B.l*B.xdd_max;
% Z_c.Mx_min1 = B.m*B.l*B.ydd_min; Z_c.Mx_max1 = -B.m*B.l*B.ydd_max;
% Z_c.Mx_min2 = B.m*B.l*B.ydd_min; Z_c.Mx_max2 = -B.m*B.l*B.ydd_max;
% Z_c.Mx_min3 = B.m*B.l*B.ydd_min; Z_c.Mx_max3 = -B.m*B.l*B.ydd_max;
% Z_c.Pe_min = D.eta*U_c.Tg_min*Z_c.omega_r_min;
% Z_c.Pe_max = -D.eta*U_c.Tg_max*Z_c.omega_r_max;
% Z_c.vr_min = 4;
% Z_c.vr_max = -25;
% Z_c.azim_min = 0;
% Z_c.azim_max = -0;

Zk = size(struct2table(Z_c),2)/2;

