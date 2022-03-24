clear all
close all

H = 144.582+4.35; % Hub height
r_base = 10/2; % Hub base outer radius
r_top = 6.5/2; % Hub top outer radius
v_m = 10; % 5, 10, 15, 20m/s
x_h = 10.93; % Hub overhang
l_b = 117.18; % Blade length
alpha = 0.15; % Wind shear exponent
R = 120.998; % Rotor radius

r_t = ((r_top-r_base)*(H-l_b))/H + r_base;

i=1;

for r_r=10:20:120
for theta=90:0.1:270
    v_i = v_m*(r_t^2*(((r_r*sind(theta))^2-x_h^2)/(x_h^2+(r_r*sind(theta))^2)^2)+((r_r*cosd(theta)+H)/H)^alpha);
    v(i) = v_i;
    i = i+1;
end
i=1;
figure(1)
plot(90:0.1:270,v)
hold on
end
xlabel('Azimuthal angle (degrees)')
ylabel('Wind speed in front of the tower (with v_m=10m/s)')
legend('R=10m','R=30m','R=50m','R=70m','R=90m','R=110m')
hold off
