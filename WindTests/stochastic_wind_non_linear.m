clear all
close all

N = 12000;  
v_m = zeros(N,1);
v_t = zeros(N,1);
Ts = 0.05; % Sampling time
ti = 0.15; % Turbulence intensity
q = 2^2/600; % Incremental variance mean wind speed
mu_m = 6; % Fixed mean wind speed: 10 m/s
L = 340.2;
sigma_m = sqrt(q); % Standard deviation mean wind noise


v_m(1) = mu_m;
v_t(1) = 0;
for i = 2:N
  w_p = (v_m(i-1)*pi)/(2*L); % Kaimal spectrum peak freq.
  %a = exp(-w_p*Ts); % Discretizing filter with zoh
  a = 1-w_p*Ts; %Discretizing filter with Fordward Euler
  sigma_t = ti*v_m(i-1)*sqrt((1-a^2)/(1-a)^2); % Standard deviation turbulent wind noise
  %sigma_t = ti*mu_m*sqrt(1-a^2);
  v_m(i) = v_m(i-1)+Ts*normrnd(0,sigma_m);
  v_t(i) = a*v_t(i-1)+(1-a)*normrnd(0,sigma_t);
end
plot(v_m)
hold on
plot(v_t)

v_e = v_m + v_t;
plot(v_e)
xlabel('iterations')
ylabel('wind speed (m/s)')
legend('v_m','v_t','v_e')
hold off