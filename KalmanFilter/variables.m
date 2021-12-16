%% Drive train model constants
Jr = 321699000; % Rotor moment of inertia
Jg = 3.223e6; % Generator moment of inertia
% cd = 0.005; % Drive train damping
% kd = 1.409e10; % Drive train torsion stiffness
mu_d = 0.05; % Drive train mechanical losses (friction)
eta_g = 0.93; % Generator efficiency

%% Tower model constants
mt = 2475680; % Tower mass
ct = 0.005; % Tower damping
kt = 1086002; % Tower stiffness

%% Aerodynamic model constants
rho = 1.225; % Density of the air
r = 241.996/2; % Rotor radius
Ar = pi*r^2; % Rotor area

%% Wind model constants
Ts = 0.05; % Sampling time
ti = 0.1; % Turbulence intensity
q = 2^2/600; % Incremental variance mean wind speed
mu_m = 6; % Fixed mean wind speed: 10 m/s
L = 340.2;
w_p = (mu_m*pi)/(2*L); % Kaimal spectrum peak freq.
%a = exp(-w_p*Ts); % Discretizing filter with zoh
a = 1-w_p*Ts; %Discretizing filter with Fordward Euler
sigma_m = sqrt(Ts*q); % Standard deviation mean wind noise
sigma_t = ti*mu_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind noise

%% Load measured data
data1 = load('..\Bladed\DLC12_06p0_Y000_S0201').DLC12_06p0_Y000_S0201;
load('..\Bladed\performancemap_data.mat')
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
Pe = data1.Data(:,20); % Electrical power
beta = data1.Data(:,30); % Mean pitch angle (collective pitch)

u = [beta Pe]';

%% Measurements
omega_g = data1.Data(:,13); % Generator speed
omega_r = data1.Data(:,10); % Rotor speed
vr = data1.Data(:,26); % Wind speed magnitud at the hub
yt_ddot = data1.Data(:,237); % Tower fore-aft acceleration

y_me = [omega_r vr yt_ddot]';

%% Initial state vector
yt_dot1 = data1.Data(1,231);
yt1 = data1.Data(1,225);
vm1 = data1.Data(1,59);

x_i = [omega_r(1) yt_dot1 yt1 0 vm1];

x_vl = {'$\omega_r$', '$\dot{y}_t$', '$y_t$', '$v_t$', '$v_m$'};
y_vl = {'$\omega_r$', '$v_r$', '$\ddot{y}_t$'};
x_ul = {'$Angular\, velocity [\frac{rad}{s}]$', '$Velocity [\frac{m}{s}]$', ...
        '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
        '$Velocity [\frac{m}{s}]$'};


