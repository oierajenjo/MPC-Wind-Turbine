clear all
close all

H = 150; % Hub height
r_base = 10/2; % Hub base outer radius
r_top = 6.5/2; % Hub top outer radius
v_m = 10; % 5, 10, 15, 20m/s
x = 11.35; % Hub overhang
% r_r = 120; % Blade length

h0 = 100;
r = ((r_top-r_base)*h0)/H + r_base;

i=1;

for r_r=10:20:120
for theta=90:0.1:270
    v_i = v_m*r^2*(((r_r*sind(theta))^2-x^2)/(x^2+(r_r*sind(theta))^2)^2)+v_m;
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
legend('r_r=10m','r_r=30m','r_r=50m','r_r=70m','r_r=90m','r_r=110m')
hold off
