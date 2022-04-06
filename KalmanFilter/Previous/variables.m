%% Drive train model constants
Jr = 321699000; % Rotor moment of inertia
Jg = 3.223e6; % Generator moment of inertia
% cd = 0.005; % Drive train damping
% kd = 1.409e10; % Drive train torsion stiffness
mu_d = 0.05; % Drive train mechanical losses (friction)
eta_g = 0.93; % Generator efficiency

%% Tower model constants
mn = 630888; % Nacelle mass
mt = 1086002; % Tower mass
mr = 387198; % Rotor mass
m = mn + mt/3; % Tower equivalent mass
% To.m = 2475680-B.m*B.B; % Tower mass
dt = 0.005; % Tower damping ratio
ft = 0.18; % Tower freq. flapwise
ct = dt*2*mt*2*pi*ft; % Tower damping
kt = (2*pi*ft)^2*mt; % Tower stiffness
h = 144.582; % Tower height
H = h + 4.34799; % Hub height

%% Aerodynamic model constants
rho = 1.225; % Density of the air
r = 241.996/2; % Rotor radius
Ar = pi*r^2; % Rotor area

%% Wind model constants
Ts = 0.05; % Sampling time
ti = 0.15; % Turbulence intensity
q = 2^2/600; % Incremental variance mean wind speed
mu_m = 6; % Fixed mean wind speed: 10 m/s
L = 340.2;
% w_p = (mu_m*pi)/(2*L); % Kaimal spectrum peak freq.
%a = exp(-w_p*Ts); % Discretizing filter with zoh
% a = 1-w_p*Ts; %Discretizing filter with Fordward Euler
% sigma_m = sqrt(Ts*q); % Standard deviation mean wind noise
% sigma_t = ti*mu_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind noise

%% Actuator constants
Ac.omega = 2.4*pi; % Natural frequency of pitch actuator model
Ac.xi = 0.8; % Damping factor of pitch actuator model
Ac.tau = 0.1; % Generator time constant

%% Measurement constants
M.sigma_enc = 0.017;
M.sigma_acc = 0.04;
M.sigma_root = 0.01;
M.sigma_pow = 0.035;
M.sigma_vane = 1;
M.sigma_azim = 0.01;

%% Load measured data
data = load('Bladed\DLC12_06p0_Y000_S0201').DLC12_06p0_Y000_S0201;
load('Bladed\performancemap_data.mat')
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
% Pe = data1.Data(:,28); % Electrical power
% beta = data1.Data(:,30); % Mean pitch angle (collective pitch)
% 
% u = [beta Pe]';
theta_ref = data.Data(:,30); % Mean pitch angle (collective pitch)
tg_ref = data.Data(:,20); % Generator Torque

u = [theta_ref tg_ref]';

%% Disturbances
vm = data.Data(:,59); % Wind mean speed

d = vm';

%% Measurements
omega_r = data.Data(:,10); % Rotor speed
xt_ddot = data.Data(:,236); % Tower fore-aft acceleration
yt_ddot = data.Data(:,237); % Tower edgewise acceleration
% Mx = sum([data1.Data(:,111) data1.Data(:,119) data1.Data(:,127)], 2);
My = sum([data1.Data(:,112) data1.Data(:,120) data1.Data(:,128)], 2);
Mx = sum([data.Data(:,61) data.Data(:,69) data.Data(:,77)], 2); % Mx in the principal axis
% My = sum([data.Data(:,62) data.Data(:,70) data.Data(:,78)], 2); % My in the principal axis
Pe = data.Data(:,28);
vr = data.Data(:,26); % Wind speed magnitud at the hub
psi = data.Data(:,11);

y_me = [omega_r vr xt_ddot]';

%% Initial state vector
xt_dot = data.Data(1,230);
xt = data.Data(1,224);
yt_dot = data.Data(1,231);
yt = data.Data(1,225);
xb = mean([data.Data(1,85) data.Data(1,91) data.Data(1,97)], 2);
xb_dot = 0;
yb = mean([data.Data(1,86) data.Data(1,92) data.Data(1,98)], 2);
yb_dot = 0;
theta = theta_ref(1);
theta_dot = mean(data.Data(1,37:39), 2);
Tg = tg_ref(1);
vt = 0;
% vm = data1.Data(1,59);

x_i = [omega_r(1) xt xt_dot yt yt_dot theta theta_dot Tg vt];

x_vl = {'$\omega_r$', '$\dot{y}_t$', '$y_t$', '$v_t$', '$v_m$'};
y_vl = {'$\omega_r$', '$v_r$', '$\ddot{y}_t$'};
x_ul = {'$Angular\, velocity [\frac{rad}{s}]$', '$Velocity [\frac{m}{s}]$', ...
        '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
        '$Velocity [\frac{m}{s}]$'};


