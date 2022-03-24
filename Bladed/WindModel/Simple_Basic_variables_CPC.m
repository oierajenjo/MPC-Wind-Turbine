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
tg_ref = data.Data(:,20); % Generator Torque
u_b = tg_ref';

%% Measurements
yt_ddot = data.Data(:,237); % Tower edgewise acceleration
y_me = yt_ddot';

%% Initial state vector
yt = data.Data(1,225);
yt_dot = data.Data(1,231);
Tg = tg_ref(1);
x_i = [yt_dot yt Tg]';

clearvars -except x_i y_me u_b d_b N data

%% Plotting variables
x_vl = {'$\dot{y}_t$', '$y_t$', '$T_g$'};

y_vl = {'$\ddot{y}_t$'};

x_ul = {'$Velocity [\frac{m}{s}]$', '$Position [m]$','$Torque [Nm]$'};