%% Drive train model constants
D.Jr = 321699000; % Rotor moment of inertia
D.Jg = 3.223e6; % Generator moment of inertia
% D.c = 0.005; % Drive train damping
% D.k = 1.409e10; % Drive train torsion stiffness
D.mu = 0; % Drive train mechanical losses (friction)
D.eta = 0.93; % Generator efficiency

%% Blades model constants
B.m = 65566/3; % Blade mass
B.d = 0.03; % Blade damping ratio
B.fx = 0.541; % Blade freq. flapwise
B.fy = 0.636; % Blade freq. edgewise
B.cx = B.d*2*B.m*2*pi*B.fx; % Blade damping x direction
B.kx = (2*pi*B.fx)^2*B.m; % Blade stiffness x direction
B.cy = B.d*2*B.m*2*pi*B.fy; % Blade damping y direction
B.ky = (2*pi*B.fy)^2*B.m; % Blade stiffness y direction
B.l = 117.1836; % Blade length
B.B = 3; % Blade amount

%% Tower model constants
Mn = 630888; % Nacelle mass
Mt = 1086002; % Tower mass
To.m = Mn + Mt/3; % Tower mass
% To.m = 2475680-B.m*B.B; % Tower mass
To.d = 0.005; % Tower damping ratio
To.f = 0.18; % Tower freq. flapwise
To.c = To.d*2*To.m*2*pi*To.f; % Tower damping
To.k = (2*pi*To.f)^2*To.m; % Tower stiffness
To.h = 144.582; % Tower height
To.r_top = 3.25; % Tower top radius
To.r_base = 5; % Tower base radius
To.H = To.h + 4.34799; % Hub height
To.r = (To.r_top-To.r_base)*(To.H-B.l)/To.H + To.r_base; % Tower radius
To.xh = 10.93; % Hub overhang

%% Aerodynamic model constants
Ae.rho = 1.225; % Density of the air
Ae.Rr = 241.996/2; % Rotor radius
Ae.Ar = pi*Ae.Rr^2; % Rotor area

%% Wind model constants
Ts = 0.05; % Sampling time
W.ti = 0.1; % Turbulence intensity
W.q = 2^2/600; % Incremental variance mean wind speed
W.mu_m = 6; % Fixed mean wind speed: 10 m/s
W.L = 340.2;
W.alpha = 0.15; % Wind shear exponent for smooth terrain

%% Actuator constants
Ac.omega = 2.4*pi; % Natural frequency of pitch actuator model
Ac.xi = 0.8; % Damping factor of pitch actuator model
Ac.tau = 0.1; % Generator time constant

%% Measurement constants
M.sigma_enc = 0.01;
M.sigma_acc = 0.02;
M.sigma_root = 0.01; % ¿?
M.sigma_pow = 0.01; % ¿?
M.sigma_vane = 1;
M.sigma_azim = 0.01;

%% Load measured data
data = load('BladedFiles\DLC12_06p0_Y000_S0201').DLC12_06p0_Y000_S0201;
% data = load('BladedFiles\DLC12_08p0_Y000_S0301').DLC12_08p0_Y000_S0301;
% data = load('BladedFiles\DLC12_10p0_Y000_S0401').DLC12_10p0_Y000_S0401;
% data = load('BladedFiles\DLC12_12p0_Y000_S0501').DLC12_12p0_Y000_S0501;
% data = load('BladedFiles\DLC12_14p0_Y000_S0601').DLC12_14p0_Y000_S0601;
% data = load('BladedFiles\DLC12_16p0_Y000_S0701').DLC12_16p0_Y000_S0701;
% data = load('BladedFiles\DLC12_18p0_Y000_S0801').DLC12_18p0_Y000_S0801;
% data = load('BladedFiles\DLC12_20p0_Y000_S0901').DLC12_20p0_Y000_S0901;
% data = load('BladedFiles\DLC12_22p0_Y000_S1001').DLC12_22p0_Y000_S1001;
% data = load('BladedFiles\DLC12_24p0_Y000_S1101').DLC12_24p0_Y000_S1101;

N = data.Channels.Scans; % Number of time steps for filter

%% Inputs
theta_ref = data.Data(:,30); % Mean pitch angle (collective pitch)
tg_ref = data.Data(:,20); % Generator Torque

u_b = [theta_ref tg_ref]';

%% Measurements
omega_r = data.Data(:,10); % Rotor speed
xt_ddot = -data.Data(:,236); % Tower fore-aft acceleration
yt_ddot = data.Data(:,237); % Tower edgewise acceleration
Pe = data.Data(:,28);
vr = data.Data(:,26); % Wind speed magnitud at the hub

y_me = [omega_r xt_ddot yt_ddot Pe vr]';

%% Initial state vector
xt_dot = -data.Data(1,230);
xt = -data.Data(1,224);
yt_dot = data.Data(1,231);
yt = data.Data(1,225);

theta = theta_ref(1);
theta_dot = mean(data.Data(1,37:39), 2);
Tg = tg_ref(1);
vt = 0;
vm = data.Data(1,59);

x_i = [omega_r(1) xt xt_dot yt yt_dot theta theta_dot Tg vt vm]';

clearvars -except D To B Ae Ac M Ts W w_p x_i y_me u_b N data

load('BladedFiles\performancemap_data.mat')
%% Plotting variables
x_vl = {'$\omega_r$', '$\dot{x}_t$', '$x_t$', '$\dot{y}_t$', '$y_t$', ...
    '$\dot{x}_{b}$', '$x_{b}$', '$\dot{y}_{b}$', '$y_{b}$',...
    '$\dot{\theta}_{b}$', '$\theta_{b}$', '$T_g$', '$v_t$', '$v_m$'};

y_vl = {'$\omega_r$', '$\ddot{x}_t$', '$\ddot{y}_t$', '$M_{y}$', '$M_{x}$',...
    '$P_e$', '$v_r$'};

x_ul = {'$Angular\, velocity [\frac{rad}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Angular\, velocity [\frac{rad}{s}]$', '$Angular\, position [rad]$',...
    '$Torque [Nm]$', '$Velocity [\frac{m}{s}]$', '$Velocity [\frac{m}{s}]$'};