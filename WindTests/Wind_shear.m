clear all
close all

H = 144.582+4.35; % Hub height
r_base = 10/2; % Hub base outer radius
r_top = 6.5/2; % Hub top outer radius
v_m = 10; % 5, 10, 15, 20m/s
x_h = 10.93; % Hub overhang
R = 120.998; % Rotor radius
l_b = 117.18; % Blade length
alpha = 0.15; % Wind shear exponent

i=1;
r_r = 80;
for r_r=10:20:120
for theta=0:0.1:360
    % w_i = alpha*(r_r/H)*cosd(theta)+alpha*(alpha-1)/2*(r_r/H)^2*cosd(theta)^2+alpha*(alpha-1)*(alpha-2)/6*(r_r/H)^3*cosd(theta)^3;
    % v_i = v_m*(1+w_i);
    v_i = v_m*((r_r*cosd(theta)+H)/H)^alpha;
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
legend('R=10m','R=30m','R=50m','R=70m','R=90m','R=110m')
hold off