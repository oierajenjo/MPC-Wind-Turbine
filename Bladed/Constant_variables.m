%% Drive train model constants
D.Jr = 321699000; % Rotor moment of inertia
D.Jg = 3.223e6; % Generator moment of inertia
% D.c = 0.005; % Drive train damping
% D.k = 1.409e10; % Drive train torsion stiffness
D.mu = 0.05; % Drive train mechanical losses (friction)
D.eta = 0.93; % Generator efficiency

%% Blades model constants
B.Mb = 65566; % Blade mass
B.m = B.Mb; % Blade equivalent mass
B.d = 0.03; % Blade damping ratio
B.fx = 0.541*sqrt(1-B.d^2); % Blade freq. flapwise
B.fy = 0.636*sqrt(1-B.d^2); % Blade freq. edgewise
B.fx = 0.541; % Blade freq. flapwise
B.fy = 0.636; % Blade freq. edgewise
B.wx = 2*pi*B.fx; % Blade x direction natural pulsation
B.wy = 2*pi*B.fy; % Blade y direction natural pulsation
B.cx = B.d*2*B.Mb*B.wx; % Blade damping x direction
B.kx = B.wx^2*B.Mb; % Blade stiffness x direction
B.cy = B.d*2*B.Mb*B.wy; % Blade damping y direction
B.ky = B.wy^2*B.Mb; % Blade stiffness y direction
B.l = 117.1836; % Blade length
B.B = 3; % Blade amount
B.cc = 1/5; % Correction coefficient Fx
B.xdd_min = 0;
B.xdd_max = 0;
B.ydd_min = 0;
B.ydd_max = 0;

%% Tower model constants
To.Mn = 630888; % Nacelle mass
To.Mt = 1086002; % Tower mass
To.Mr = 387198; % Rotor mass
To.m = To.Mn + To.Mt/3; % Tower equivalent mass
% To.m = 2475680-B.m*B.B; % Tower mass
To.d = 0.005; % Tower damping ratio
To.f = 0.18; % Tower freq. flapwise
To.c = To.d*2*To.Mt*2*pi*To.f; % Tower damping
To.k = (2*pi*To.f)^2*To.Mt; % Tower stiffness
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

%% Aerodynamic model constants
Ae.rho = 1.225; % Density of the air
Ae.Rr = 241.996/2; % Rotor radius
Ae.Ar = pi*Ae.Rr^2; % Rotor area

%% Wind model constants
Ts = 0.05; % Sampling time
W.ti = 0.15; % Turbulence intensity
W.q = 2^2/600; % Incremental variance mean wind speed
W.mu_m = 6; % Fixed mean wind speed: 10 m/s
W.L = 340.2;
W.alpha = 0.15; % Wind shear exponent for smooth terrain

%% Actuator constants
Ac.omega = 2.4*pi; % Natural frequency of pitch actuator model
Ac.xi = 0.8; % Damping factor of pitch actuator model
Ac.tau = 0.1; % Generator time constant

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

%% Vector Sizes
Lk = size(x_i,1); % Size of state vector
Yk = size(y_me,1); % Size of measured vector
Uk = size(u_b,1); % Size of imput vector
t = Ts*(1:N);