close all
clear all

data = load('Bladed\DLC12_06p0_Y000_S0201').DLC12_06p0_Y000_S0201;
load('Bladed\performancemap_data.mat')
ts = data.Data(:,25); % Generator Torque
theta_ref = data.Data(:,30); % Mean pitch angle (collective pitch)
tg_ref = data.Data(:,20); % Generator Torque
vm = data.Data(:,59); % Wind mean speed
omega_r = data.Data(:,10); % Rotor speed

m_b = 65566;
B = 3;
c_bx = 13372;
k_bx = 757590;
c_by = 15721;
k_by = 1047000;
m_t = 2475680-B*m_b;
c_t = 25775;
k_t = 2915000;
H = 148.93;
rho = 1.225; % Density of the air
Rr = 241.996/2; % Rotor radius
Ar = pi*Rr^2; % Rotor area

for j=1:2
    for i=1:length(ts)
        u(i,j) = 1;
    end
end
% u = u';

for i=1:length(ts)
    Fr(i) = 0.5*rho*Ar*vm(i)^2*cp_ct(omega_r(i)*Rr/(vm(i)),theta_ref(i),ct_l,lambdaVec,pitchVec);
end
Fr = Fr';

figure(20)
plot(Fr/(B*m_b))
legend('Fr/(B*m_b)')
figure(21)
plot((3*tg_ref)/(2*H*m_t))
legend('(3*tg_{ref})/(2*H*m_t)')

% Initial conditions
xt_dot_i = data.Data(1,230);
xt_i = data.Data(1,224);
xb_i = mean([data.Data(1,85) data.Data(1,91) data.Data(1,97)], 2);
xb_dot_i = 0;
yt_dot_i = data.Data(1,231);
yt_i = data.Data(1,225);
yb_i = mean([data.Data(1,86) data.Data(1,92) data.Data(1,98)], 2);
yb_dot_i = 0;

% Tower and blades flapwise
A1 = [0 1 0 0;
    (-B*k_bx-k_t)/m_t (-B*c_bx-c_t)/m_t B*k_bx/m_t B*c_bx/m_t;
    0 0 0 1;
    k_bx/m_b c_bx/m_b -k_bx/m_b -c_bx/m_b];
B1 = [0 0 0 1/(B*m_b)]';
C1 = eye(4);
D1 = zeros(4,1);

sys1=ss(A1,B1,C1,D1);

x_i_1 = [xt_i xt_dot_i xb_i xb_dot_i];

figure(1)
[y_1,ts_1,x_1] = lsim(sys1,Fr,ts,x_i_1);
plot(ts_1,y_1)
xlabel('Time (sec)')
ylabel('System response')
legend('xt','xt_{dot}','xb','xb_{dot}')

% Tower and blades edgewise
A2 = [0 1 0 0;
    (-B*k_by-k_t)/m_t (-B*c_by-c_t)/m_t B*k_by/m_t B*c_by/m_t;
    0 0 0 1;
    k_by/m_b c_by/m_b -k_by/m_b -c_by/m_b];
B2 = [0 3/(2*H*m_t) 0 0]';
C2 = eye(4);
D2 = zeros(4,1);

sys2=ss(A2,B2,C2,D2);

x_i_2 = [yt_i yt_dot_i yb_i yb_dot_i];

figure(2)
[y_2,ts_2,x_2] = lsim(sys2,tg_ref,ts,x_i_2);
plot(ts_2,y_2)
xlabel('Time (sec)')
ylabel('System response')
legend('yt','yt_{dot}','yb','yb_{dot}')

% Tower fore-aft
A3 = [0 1;
    -k_t/m_t -c_t/m_t];
B3 = [0 1/m_t]';
C3 = eye(2);
D3 = zeros(2,1);

sys3=ss(A3,B3,C3,D3);

x_i_3 = [xt_i xt_dot_i];

figure(3)
[y_3,ts_3,x_3] = lsim(sys3,Fr,ts,x_i_3);
plot(ts_3,y_3)
xlabel('Time (sec)')
ylabel('System response')
legend('xt','xt_{dot}')

% Tower sidewards
A4 = [0 1;
    -k_t/m_t -c_t/m_t];
B4 = [0 3/(2*H*m_t)]';
C4 = eye(2);
D4 = zeros(2,1);

sys4=ss(A4,B4,C4,D4);

x_i_4 = [yt_i yt_dot_i];

figure(4)
[y_4,ts_4,x_4] = lsim(sys4,tg_ref,ts,x_i_4);
plot(ts_4,y_4)
xlabel('Time (sec)')
ylabel('System response')
legend('yt','yt_{dot}')

% Tower and bldes flapwise and sideways
A5 = [0 1 0 0 0 0 0 0;
    (-B*k_bx-k_t)/m_t (-B*c_bx-c_t)/m_t 0 0 B*k_bx/m_t B*c_bx/m_t 0 0;
    0 0 0 1 0 0 0 0;
    (-B*k_by-k_t)/m_t (-B*c_by-c_t)/m_t 0 0 B*k_by/m_t B*c_by/m_t 0 0;
    0 0 0 0 0 1 0 0;
    k_bx/m_b c_bx/m_b 0 0 -k_bx/m_b -c_bx/m_b 0 0;
    0 0 0 0 0 0 0 1;
    0 0 k_by/m_b c_by/m_b 0 0 -k_by/m_b -c_by/m_b];
B5 = [0 0 0 0 0 1/(B*m_b) 0 0;
    0 0 0 3/(2*H*m_t) 0 0 0 0]';
C5 = eye(8);
D5 = zeros(8,2);

sys5=ss(A5,B5,C5,D5);

x_i_5 = [xt_i xt_dot_i yt_i yt_dot_i xb_i xb_dot_i yb_i yb_dot_i];

figure(5)
[y_5,ts_5,x_5] = lsim(sys5,u',ts,x_i_5);
plot(ts_5,y_5,'Linewidth',2)
xlabel('Time (sec)')
ylabel('System response')
legend('xt','xt_{dot}','yt','yt_{dot}','xb','xb_{dot}','yb','yb_{dot}')


function res = cp_ct(la,be,cl,lambdaVec,pitchVec)
[~,i_la] = min(abs(lambdaVec-abs(la)));
[~,i_be] = min(abs(pitchVec-be));
res = cl(i_la,i_be);
end