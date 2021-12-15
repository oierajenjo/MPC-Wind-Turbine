close all
clear all

% Drive train model constants
Jr = 321699000; % Rotor moment of inertia
Jg = 3.223e6; % Generator moment of inertia
cd = 0.005; % Drive train damping
kd = 1.409e6; % Drive train torsion stiffness
mu_d = 0.05; % Drive train mechanical losses (friction)
eta_g = 0.93; % Generator efficiency
% Tower model constants
mt = 1086002; % Tower mass
ct = 0.005; % Tower damping
kt = 2.028635714e12; % Tower stiffness (avg)
% Aerodynamic model constants
rho = 1.225; % Density of the air
R = 241.996/2; % Rotor radius
Ar = pi*R^2; % Rotor area
% Wind model constants
Ts = 0.05; % Sampling time
ti = 0.2; % Turbulence intensity
q = 2^2/600; % Incremental variance mean wind speed
mu_m = 10; % Fixed mean wind speed: 10 m/s
L = 340.2;
w_p = (mu_m*pi)/(2*L); % Kaimal spectrum peak freq.
%a = exp(-w_p*Ts); % Discretizing filter with zoh
a = 1-w_p*Ts; %Discretizing filter with Fordward Euler
sigma_m = sqrt(Ts*q); % Standard deviation mean wind noise
sigma_t = ti*mu_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind noise

% Load data
data1 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_06p0_Y000_S0201');
data2 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_08p0_Y000_S0301');
data3 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_10p0_Y000_S0401');
data4 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_12p0_Y000_S0501');
data5 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_14p0_Y000_S0601');
data6 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_16p0_Y000_S0701');
data7 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_18p0_Y000_S0801');
data8 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_20p0_Y000_S0901');
data9 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_22p0_Y000_S1001');
data10 = load('C:\Alvaro\AAU Master\SEMESTER 3-4\THESIS\Bladed data\DLC12_24p0_Y000_S1101');
% Inputs
Pe = data1.DLC12_06p0_Y000_S0201.Data(:,20); % Electrical power
beta = data1.DLC12_06p0_Y000_S0201.Data(:,30); % Mean pitch angle (collective pitch)
% Measurements
omega_g = data1.DLC12_06p0_Y000_S0201.Data(:,13); % Generator speed
omega_r = data1.DLC12_06p0_Y000_S0201.Data(:,10); % Rotor speed
vr = data1.DLC12_06p0_Y000_S0201.Data(:,26); % Wind speed magnitud at the hub
yt_ddot = data1.DLC12_06p0_Y000_S0201.Data(:,237); % Tower fore-aft acceleration
% Extra available data
% yt_dot = data1.DLC12_06p0_Y000_S0201.Data(:,231); % Tower fore-aft velocity
% yt = data1.DLC12_06p0_Y000_S0201.Data(:,225); % Tower fore-aft deflection

% Euler's Discretization Method
% Initial conditions and setup
n = 12000/Ts;
yt(1) = data1.DLC12_06p0_Y000_S0201.Data(1,225); % Tower fore-aft velocity initial value
yt_dot(1) = data1.DLC12_06p0_Y000_S0201.Data(1,231); % Tower fore-aft deflection initial value

% x = (enter the starting value of x here):Ts:(enter the ending value of x here);  % the range of x
% y = zeros(size(x));  % allocate the result y
% y(1) = (enter the starting value of y here);  % the initial y value
% n = numel(y);  % the number of y values
% The loop to solve the differential equations
for i=1:n-1
    % Wind
    nm(i+1) = normrnd(0,sigma_m^2);
    nt(i+1) = normrnd(0,sigma_t^2);
    vm(i+1) = vm(i)+nm(i);
    vt(i+1) = a*vt(i)+(1-a)*nt(i);
    
    yt(i+1) = yt(i)+Ts*(vm(i)+vt(i)-vr(i));
    
    % Aerodynamics
    lambda(i) = (omega_r(i)*R)/vr(i);
    Tr(i) = 1/2*rho*vr(i)^3*Ar*Cp(lambda(i),beta(i))*1/omega_r(i);
    Fr(i) = 1/2*rho*vr(i)^2*Ar*Ct(lambda(i),beta(i));
    Tg(i) = Pe(i)/(eta_g*omega_g(i));
    
    % Drive train (two inertia)
%     gamma_dot(i) = omega_r(i)-omega_g(i);
%     gamma(i) = trapz(gamma_dot);
%     omega_r_dot(i+1) = omega_r_dot(i)+Ts*((1-mu_d)*1/Jr*Tr-cd/Jr*gamma_dot(i)-kd/Jr*gamma(i));
%     omega_g_dot(i+1) = omega_g_dot(i)+Ts*(cd/Jg*gamma_dot(i)+kd/Jg*gamma(i)-Tg(i));

    % Drive train (one inertia)
    omega_r(i+1) = omega_r(i)+Ts*((1-mu_d)*Tr/Jr-Tg/Jr);
    
    % Tower
    yt_dot(i+1) = yt_dot(i)+Ts*(1/mt*Fr(i)-ct/mt*yt_dot(i)-kt/mt*yt(i));
end

function x = WindTurbStateFcn(x,ts)
    % Discrete-time approx. of the wind turbine model
    
    % Drive train model constants
    Jr = 321699000; % Rotor moment of inertia
    Jg = 3.223e6; % Generator moment of inertia
    cd = 0.005; % Drive train damping
    kd = 1.409e6; % Drive train torsion stiffness
    mu_d = 0.05; % Drive train mechanical losses (friction)
    eta_g = 0.93; % Generator efficiency
    % Tower model constants
    mt = 1086002; % Tower mass
    ct = 0.005; % Tower damping
    kt = 2.028635714e12; % Tower stiffness (avg)
    % Aerodynamic model constants
    rho = 1.225; % Density of the air
    R = 241.996/2; % Rotor radius
    Ar = pi*R^2; % Rotor area
    % Wind model constants
    Ts = 0.05; % Sampling time
    ti = 0.2; % Turbulence intensity
    q = 2^2/600; % Incremental variance mean wind speed
    mu_m = 10; % Fixed mean wind speed: 10 m/s
    L = 340.2;
    w_p = (mu_m*pi)/(2*L); % Kaimal spectrum peak freq.
    %a = exp(-w_p*Ts); % Discretizing filter with zoh
    a = 1-w_p*Ts; %Discretizing filter with Fordward Euler
    sigma_m = sqrt(Ts*q); % Standard deviation mean wind noise
    sigma_t = ti*mu_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind noise
    
    % Euler integration of continuous-time dynamics x'=f(x) with sample time dt
    omega_r = omega_r+Ts*((1-mu_d)*Tr/Jr-Tg/Jr);
end

function dxdt = WindTurbStateFcnContinuous(x)
    dxdt = [(1-mu_d)*Tr/Jr-Tg/Jr;
        ];
end
