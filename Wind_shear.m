clear all
close all

H = 150; % Hub height
r_base = 10; % Hub base outer radius
r_top = 6.5; % Hub top outer radius
v_m = 10; % 5, 10, 15, 20m/s
x = 11.35; % Hub overhang
% r_r = 120; % Blade length
alpha = 0.1; % Wind shear exponent

i=1;

for r_r=10:20:120
for theta=0:0.1:360
    w_i = alpha*(r_r/H)*cosd(theta)+alpha*(alpha-1)/2*(r_r/H)^2*cosd(theta)^2+alpha*(alpha-1)*(alpha-2)/6*(r_r/H)^3*cosd(theta)^3;
    v_i = v_m*(1+w_i);
    v(i) = v_i;
    i = i+1;
end
i=1;
figure(1)
plot(0:0.1:360,v)
hold on
end
xlabel('Azimuthal angle (degrees)')
ylabel('Wind speed (with v_m=10m/s)')
legend('r_r=10m','r_r=30m','r_r=50m','r_r=70m','r_r=90m','r_r=110m')
hold off