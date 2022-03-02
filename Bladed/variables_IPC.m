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

%% Tower model constants
To.m = 2475680; % Tower mass
% T.m = 1389678; % Nacelle mass
To.d = 0.005; % Tower damping ratio
To.f = 0.18; % Tower freq. flapwise
To.c = To.d*2*To.m*2*pi*To.f; 
To.k = (2*pi*To.f)^2*To.m; % Tower stiffness
To.h = 144.582; % Tower height
To.r_top = 3.25; % Tower top radius
To.r_base = 5; % Tower base radius
To.H = To.h + 4.34799; % Hub height
To.r = (To.r_top-To.r_base)*(To.H-B.l)/To.H + To.r_base; % Tower radius
To.xh = 10.93; % Tower overhang

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
data = load('DLC12_06p0_Y000_S0201').DLC12_06p0_Y000_S0201;
% data = load('..\Bladed\DLC12_08p0_Y000_S0301').DLC12_08p0_Y000_S0301;
% data = load('..\Bladed\DLC12_10p0_Y000_S0401').DLC12_10p0_Y000_S0401;
% data = load('..\Bladed\DLC12_12p0_Y000_S0501').DLC12_12p0_Y000_S0501;
% data = load('..\Bladed\DLC12_14p0_Y000_S0601').DLC12_14p0_Y000_S0601;
% data = load('..\Bladed\DLC12_16p0_Y000_S0701').DLC12_16p0_Y000_S0701;
% data = load('..\Bladed\DLC12_18p0_Y000_S0801').DLC12_18p0_Y000_S0801;
% data = load('..\Bladed\DLC12_20p0_Y000_S0901').DLC12_20p0_Y000_S0901;
% data = load('..\Bladed\DLC12_22p0_Y000_S1001').DLC12_22p0_Y000_S1001;
% data = load('..\Bladed\DLC12_24p0_Y000_S1101').DLC12_24p0_Y000_S1101;ç

% time = data.Data(:,1);

%% Inputs
theta_ref = data.Data(:,34:36); % Mean pitch angle (collective pitch)
tg_ref = data.Data(:,20); % Generator Torque

u = [theta_ref tg_ref]';

%% Measurements
omega_r = data.Data(:,10); % Rotor speed
xt_ddot = data.Data(:,236); % Tower fore-aft acceleration
yt_ddot = data.Data(:,237); % Tower edgewise acceleration
Mx = [data.Data(:,111) data.Data(:,119) data.Data(:,127)];
My = [data.Data(:,112) data.Data(:,120) data.Data(:,128)];
% Mx = [data.Data(:,61) data.Data(:,69) data.Data(:,77)]; % Mx in the principal axis
% My = [data.Data(:,62) data.Data(:,70) data.Data(:,78)]; % My in the principal axis
Pe = data.Data(:,28);
vr = data.Data(:,26); % Wind speed magnitud at the hub
psi = data.Data(:,11);

y_me = [omega_r xt_ddot yt_ddot My Mx Pe vr psi]';

%% Initial state vector
xt_dot = data.Data(1,230);
xt = data.Data(1,224);
yt_dot = data.Data(1,231);
yt = data.Data(1,225);
xb = [data.Data(1,85) data.Data(1,91) data.Data(1,97)];
yb = [data.Data(1,86) data.Data(1,92) data.Data(1,98)];
xb_dot = zeros(1,3);
yb_dot = zeros(1,3);
theta = theta_ref(1,:);
theta_dot = data.Data(1,37:39);
Tg = tg_ref(1);
vt = 0;
vm = data.Data(1,59);

x_i = [omega_r(1) xt xt_dot yt yt_dot xb xb_dot yb yb_dot theta theta_dot Tg vt vm psi(1)]';

N = data.Channels.Scans; % Number of time steps for filter
clearvars -except D To B Ae Ac M Ts ti W w_p x_i y_me u N % a sigma_m sigma_t

load('performancemap_data.mat')
%% Plotting variables
x_vl = {'$\omega_r$', '$\dot{x}_t$', '$x_t$', '$\dot{y}_t$', '$y_t$', ...
    '$\dot{x}_{b_1}$', '$\dot{x}_{b_2}$', '$\dot{x}_{b_3}$',...
    '$x__{b_1}$', '$x__{b_2}$', '$x__{b_3}$', '$\dot{y}_{b_1}$',...
    '$\dot{y}_{b_2}$', '$\dot{y}_{b_3}$', '$y__{b_1}$', '$y__{b_2}$',...
    '$y__{b_3}$', '$\dot{\theta}_{b_1}$', '$\dot{\theta}_{b_2}$',...
    '$\dot{\theta}_{b_3}$', '$\theta__{b_1}$', '$\theta__{b_2}$',...
    '$\theta__{b_3}$', '$T_g$', '$v_t$', '$v_m$', '$\psi$'};


y_vl = {'$\omega_r$', '$\ddot{x}_t$', '$\ddot{y}_t$', '$M_{x_1}$', '$M_{x_2}$',...
    '$M_{x_3}$', '$M_{y_1}$', '$M_{y_2}$', '$M_{y_3}$', '$P_e$', '$v_r$', '$psi$'};

x_ul = {'$Angular\, velocity [\frac{rad}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Velocity [\frac{m}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Position [m]$', '$Position [m]$', ...
    '$Velocity [\frac{m}{s}]$', '$Velocity [\frac{m}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', ...
    '$Position [m]$', '$Position [m]$', ...
    '$Angular\, velocity [\frac{rad}{s}]$', '$Angular\, velocity [\frac{rad}{s}]$',...
    '$Angular\, velocity [\frac{rad}{s}]$', '$Angular\, position [rad]$', ...
    '$Angular\, position [rad]$', '$Angular\, position [rad]$',... 
    '$Torque [Nm]$', '$Velocity [\frac{m}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Angular\, position [rad]$'};