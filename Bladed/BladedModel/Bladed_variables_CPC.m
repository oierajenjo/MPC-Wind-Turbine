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
Ts = 0.05;

%% Inputs
theta_ref = data.Data(:,30); % Mean pitch angle (collective pitch)
tg_ref = data.Data(:,20); % Generator Torque

u_b = [theta_ref tg_ref]';

%% Disturbances
vr = data.Data(:,54); % Wind speed
d_b = vr';

%% Measurements
omega_r = data.Data(:,10); % Rotor speed
xt_ddot = data.Data(:,236); % Tower fore-aft acceleration
yt_ddot = data.Data(:,237); % Tower edgewise acceleration
% Mx = sum([data.Data(:,111) data.Data(:,119) data.Data(:,127)], 2); % Mx in the principal axis (hub)
% My = sum([data.Data(:,112) data.Data(:,120) data.Data(:,128)], 2); % My in the principal axis (hub)
Mx = sum([data.Data(:,61) data.Data(:,69) data.Data(:,77)], 2); % Mx in the principal axis (blade)
My = sum([data.Data(:,62) data.Data(:,70) data.Data(:,78)], 2); % My in the principal axis (blade)
Pe = data.Data(:,28);
vr = data.Data(:,54); % Wind speed magnitud at the hub

y_me = [omega_r xt_ddot yt_ddot My Mx Pe vr]';

%% Initial state vector
xt_dot = data.Data(1,230);
xt = data.Data(1,224);
yt_dot = data.Data(1,231);
yt = data.Data(1,225);

xb = mean([data.Data(1,85) data.Data(1,91) data.Data(1,97)], 2);
xb_dot = mean([data.Data(2,85)-data.Data(1,85) data.Data(2,91)-data.Data(1,91) data.Data(2,97)-data.Data(1,97)]/Ts, 2);
yb = mean([data.Data(1,86) data.Data(1,92) data.Data(1,98)], 2);
yb_dot = mean([data.Data(2,86)-data.Data(1,86) data.Data(2,92)-data.Data(1,92) data.Data(2,98)-data.Data(1,98)]/Ts, 2);

theta = theta_ref(1);
theta_dot = mean(data.Data(1,37:39), 2);
Tg = tg_ref(1);

x_i = [omega_r(1) xt xt_dot yt yt_dot xb xb_dot yb yb_dot theta theta_dot Tg]';

clearvars -except x_i y_me u_b d_b N data

%% Plotting variables
x_vl = {'$\omega_r$', '$x_t$', '$\dot{x}_t$', '$y_t$', '$\dot{y}_t$', ...
    '$x_{b}$', '$\dot{x}_{b}$', '$y_{b}$', '$\dot{y}_{b}$',...
    '$\theta_{b}$', '$\dot{\theta}_{b}$', '$T_g$'};

y_vl = {'$\omega_r$', '$\ddot{x}_t$', '$\ddot{y}_t$', '$M_{y}$', '$M_{x}$',...
    '$P_e$', '$v_r$'};

x_ul = {'$Angular\, velocity [\frac{rad}{s}]$', ...
    '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
    '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
    '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
    '$Position [m]$', '$Velocity [\frac{m}{s}]$', ...
    '$Angular\, position [rad]$', '$Angular\, velocity [\frac{rad}{s}]$',...
    '$Torque [Nm]$'};