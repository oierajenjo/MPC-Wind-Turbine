clear all
close all

T_end = 600;
Ts = 0.05; % Sampling time
N = T_end/Ts;
v_m = zeros(N,1);
v_t = zeros(N,1);
ti = 0.15; % Turbulence intensity
q = 2^2/600; % Incremental variance mean wind speed
mu_m = 6; % Fixed mean wind speed: 10 m/s
L = 340.2;
w_p = (mu_m*pi)/(2*L); % Kaimal spectrum peak freq.
% a = exp(-w_p*Ts); % Discretizing filter with zoh
a = 1-w_p*Ts; %Discretizing filter with Fordward Euler
sigma_m = sqrt(q); % Standard deviation mean wind noise
sigma_t = ti*mu_m*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind noise
% sigma_t = ti*mu_m*sqrt(1-a^2);

v_m(1) = mu_m;
for i = 2:N
  v_m(i) = v_m(i-1)+Ts*normrnd(0,sigma_m); 
end
plot(v_m)
hold on

v_t(1) = 0;
for i = 2:N
  v_t(i) = a*v_t(i-1)+(1-a)*normrnd(0,sigma_t);
end
plot(v_t)

v_e = v_m + v_t;
plot(v_e)
xlabel('time (sec.)')
ylabel('wind speed (m/s)')
legend('v_m','v_t','v_e')
hold off