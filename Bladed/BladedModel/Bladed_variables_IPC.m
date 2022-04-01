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
theta_ref = data.Data(:,34:36); % Mean pitch angle (collective pitch)
tg_ref = data.Data(:,20); % Generator Torque

u_b = [theta_ref tg_ref]';

%% Disturbances
vr = data.Data(:,54); % Wind speed
d_b = vr';

%% Measurements
omega_r = data.Data(:,10); % Rotor speed
xt_ddot = data.Data(:,236); % Tower fore-aft acceleration
yt_ddot = data.Data(:,237); % Tower edgewise acceleration
My = [data.Data(:,112) data.Data(:,120) data.Data(:,128)]; % My in the root axis
% Mx = [data.Data(:,111) data.Data(:,119) data.Data(:,127)]; % Mx in the root axis
% My = [data.Data(:,62) data.Data(:,70) data.Data(:,78)]; % My in the principal axis
Mx = [data.Data(:,61) data.Data(:,69) data.Data(:,77)]; % Mx in the principal axis
Pe = data.Data(:,28);
vr = data.Data(:,54); % Wind speed magnitud at the hub
psi = data.Data(:,11);

y_me = [omega_r xt_ddot yt_ddot My Mx Pe vr psi]';

%% Initial state vector
xt_dot = data.Data(1,230);
xt = data.Data(1,224);
yt_dot = data.Data(1,231);
yt = data.Data(1,225);

Ts = 0.05;
xb = [data.Data(1,85) data.Data(1,91) data.Data(1,97)];
xb_dot = [(data.Data(2,85)-data.Data(1,85))/Ts (data.Data(2,91)-data.Data(1,91))/Ts (data.Data(2,97)-data.Data(1,97))/Ts];
yb = [data.Data(1,86) data.Data(1,92) data.Data(1,98)];
yb_dot = [(data.Data(2,86)-data.Data(1,86))/Ts (data.Data(2,92)-data.Data(1,92))/Ts (data.Data(2,98)-data.Data(1,98))/Ts];

theta = theta_ref(1,:);
theta_dot = data.Data(1,37:39);
tg = tg_ref(1);
vt = 0;
vm = mean(data.Data(:,59));

x_i = [omega_r(1) xt xt_dot yt yt_dot xb xb_dot yb yb_dot theta theta_dot...
    tg psi(1)]';

%% Initial state vector
xt_dot = data.Data(:,230);
xt = data.Data(:,224);
yt_dot = data.Data(:,231);
yt = data.Data(:,225);

xb = [data.Data(:,85) data.Data(:,91) data.Data(:,97)];
xb_dot = zeros(N,1);
yb = [data.Data(:,86) data.Data(:,92) data.Data(:,98)];
yb_dot = zeros(N,1);

theta = theta_ref;
theta_dot = data.Data(:,37:39);
tg = tg_ref;

x_me = [omega_r xt xt_dot yt yt_dot xb xb_dot yb yb_dot theta theta_dot...
    tg psi(1)]';

clearvars -except x_i y_me u_b d_b N data

%% Plotting variables
x_vl = {'$\omega_r$', '$x_t$', '$\dot{x}_t$', '$y_t$', '$\dot{y}_t$', ...
    '$x_{b_1}$', '$x_{b_2}$', '$x_{b_3}$', '$\dot{x}_{b_1}$',...
    '$\dot{x}_{b_2}$', '$\dot{x}_{b_3}$','$y_{b_1}$', '$y_{b_2}$',...
    '$y_{b_3}$', '$\dot{y}_{b_1}$', '$\dot{y}_{b_2}$', '$\dot{y}_{b_3}$',...
    '$\theta_{b_1}$', '$\theta_{b_2}$', '$\theta_{b_3}$',...
    '$\dot{\theta}_{b_1}$', '$\dot{\theta}_{b_2}$',...
    '$\dot{\theta}_{b_3}$', '$T_g$', '$\psi$'};


y_vl = {'$\omega_r$', '$\ddot{x}_t$', '$\ddot{y}_t$', '$M_{x_1}$', '$M_{x_2}$',...
    '$M_{x_3}$', '$M_{y_1}$', '$M_{y_2}$', '$M_{y_3}$', '$P_e$', '$v_r$', '$psi$'};

x_ul = {'$Angular\, velocity [\frac{rad}{s}]$', ...
    '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
    '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
    '$Position [m]$', '$Position [m]$', '$Position [m]$',...
    '$Velocity [\frac{m}{s}]$', '$Velocity [\frac{m}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Position [m]$', '$Position [m]$', ...
    '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
    '$Velocity [\frac{m}{s}]$', '$Velocity [\frac{m}{s}]$', ...
    '$Angular\, position [rad]$', '$Angular\, position [rad]$', ...
    '$Angular\, position [rad]$', '$Angular\, velocity [\frac{rad}{s}]$',... 
    '$Angular\, velocity [\frac{rad}{s}]$',... 
    '$Angular\, velocity [\frac{rad}{s}]$',...
    '$Torque [Nm]$', '$Angular\, position [rad]$'};