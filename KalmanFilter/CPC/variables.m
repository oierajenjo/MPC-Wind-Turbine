%% Drive train model constants
D.Jr = 321699000; % Rotor moment of inertia
D.Jg = 3.223e6; % Generator moment of inertia
% D.c = 0.005; % Drive train damping
% D.k = 1.409e10; % Drive train torsion stiffness
D.mu = 0.05; % Drive train mechanical losses (friction)
D.eta = 0.93; % Generator efficiency

%% Blades model constants
B.m = 65566; % Blade mass
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
T.m = 2475680-B.m*B.B; % Tower mass
T.d = 0.005; % Tower damping ratio
T.f = 0.18; % Tower freq. flapwise
T.c = T.d*2*T.m*2*pi*T.f; % Tower damping
T.k = (2*pi*T.f)^2*T.m; % Tower stiffness
T.h = 144.582; % Tower height
T.r_top = 3.25; % Tower top radius
T.r_base = 5; % Tower base radius
T.H = T.h + 4.34799; % Hub height
T.r = (T.r_top-T.r_base)*(T.H-B.l)/T.H + T.r_base; % Tower radius
T.xh = 10.93; % Hub overhang

%% Aerodynamic model constants
Ae.rho = 1.225; % Density of the air
Ae.Rr = 241.996/2; % Rotor radius
Ae.Ar = pi*Ae.Rr^2; % Rotor area

%% Wind model constants
Ts = 0.05; % Sampling time
ti = 0.1; % Turbulence intensity
W.q = 2^2/600; % Incremental variance mean wind speed
W.mu_m = 6; % Fixed mean wind speed: 10 m/s
W.L = 340.2;
W.alpha = 0.15; % Wind shear exponent for smooth terrain
% w_p = (W.mu_m*pi)/(2*W.L); % Kaimal spectrum peak freq.
% W.a = exp(-W.w_p*Ts); % Discretizing filter with zoh
% a = 1-w_p*Ts; %Discretizing filter with Fordward Euler
% sigma_m = sqrt(Ts*W.q); % Standard deviation mean wind noise
% sigma_t = ti*W.mu_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind noise

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
data1 = load('DLC12_06p0_Y000_S0201').DLC12_06p0_Y000_S0201;
% data2 = load('..\Bladed\DLC12_08p0_Y000_S0301').DLC12_08p0_Y000_S0301;
% data3 = load('..\Bladed\DLC12_10p0_Y000_S0401').DLC12_10p0_Y000_S0401;
% data4 = load('..\Bladed\DLC12_12p0_Y000_S0501').DLC12_12p0_Y000_S0501;
% data5 = load('..\Bladed\DLC12_14p0_Y000_S0601').DLC12_14p0_Y000_S0601;
% data6 = load('..\Bladed\DLC12_16p0_Y000_S0701').DLC12_16p0_Y000_S0701;
% data7 = load('..\Bladed\DLC12_18p0_Y000_S0801').DLC12_18p0_Y000_S0801;
% data8 = load('..\Bladed\DLC12_20p0_Y000_S0901').DLC12_20p0_Y000_S0901;
% data9 = load('..\Bladed\DLC12_22p0_Y000_S1001').DLC12_22p0_Y000_S1001;
% data10 = load('..\Bladed\DLC12_24p0_Y000_S1101').DLC12_24p0_Y000_S1101;ç

% time = data1.Data(:,1);

%% Inputs
theta_ref = data1.Data(:,30); % Mean pitch angle (collective pitch)
tg_ref = data1.Data(:,20); % Generator Torque

u = [theta_ref tg_ref]';

%% Measurements
omega_r = data1.Data(:,10); % Rotor speed
xt_ddot = data1.Data(:,236); % Tower fore-aft acceleration
yt_ddot = data1.Data(:,237); % Tower edgewise acceleration
% Mx = mean([data1.Data(:,111) data1.Data(:,119) data1.Data(:,127)], 2);
% My = mean([data1.Data(:,112) data1.Data(:,120) data1.Data(:,128)], 2);
Mx = mean([data1.Data(:,61) data1.Data(:,69) data1.Data(:,77)], 2); % Mx in the principal axis
My = mean([data1.Data(:,62) data1.Data(:,70) data1.Data(:,78)], 2); % My in the principal axis
Pe = data1.Data(:,28);
vr = data1.Data(:,26); % Wind speed magnitud at the hub
psi = data1.Data(:,11);

y_me = [omega_r xt_ddot yt_ddot My Mx Pe vr]';

%% Initial state vector
xt_dot = data1.Data(1,230);
xt = data1.Data(1,224);
yt_dot = data1.Data(1,231);
yt = data1.Data(1,225);
xb = mean([data1.Data(1,85) data1.Data(1,91) data1.Data(1,97)], 2);
yb = mean([data1.Data(1,86) data1.Data(1,92) data1.Data(1,98)], 2);
xb_dot = 0;
yb_dot = 0;
theta = theta_ref(1);
theta_dot = mean(data1.Data(1,37:39), 2);
Tg = tg_ref(1);
vt = 0;
vm = data1.Data(1,59);

x_i = [omega_r(1) xt xt_dot yt yt_dot xb xb_dot yb yb_dot theta theta_dot Tg vt vm];

N = data1.Channels.Scans; % Number of time steps for filter
clearvars -except D T B Ae Ac M Ts ti W w_p x_i y_me u N data1% a sigma_m sigma_t

load('performancemap_data.mat')
%% Plotting variables
x_vl = {'$\omega_r$', '$\dot{x}_t$', '$x_t$', '$\dot{y}_t$', '$y_t$', ...
    '$\dot{x}_{b}$', '$x__{b}$', '$\dot{y}_{b}$', '$y__{b}$',...
    '$\dot{\theta}_{b}$', '$\theta__{b}$', '$T_g$', '$v_t$', '$v_m$', '$\psi$'};

% = [wr xt x˙t yt y˙t xb x˙b yb y˙b q ˙q Tg vt vm y]T
y_vl = {'$\omega_r$', '$\ddot{x}_t$', '$\ddot{y}_t$', '$M_{y}$', '$M_{x}$',...
    '$P_e$', '$v_r$', '$psi$'};

x_ul = {'$Angular\, velocity [\frac{rad}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Angular\, velocity [\frac{rad}{s}]$', '$Angular\, position [rad]$',...
    '$Torque [Nm]$', '$Velocity [\frac{m}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Angular\, position [rad]$'};