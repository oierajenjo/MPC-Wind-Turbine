function [] = vri_plot(var,xk,xm,t)
[~,Ae,~,W,~,~,To] = extract_struct(var);

ws_ts = @(x,i) (To.r^2*(Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2-To.xh^2)/(To.xh^2+Ae.Rr^2*(sin(x(27)+2*pi*i/3))^2)^2 +...
    ((Ae.Rr*cos(x(27)+2*pi*i/3)+To.H)/To.H)^W.alpha); % Wind Share and Tower Shadow
vei = @(x,i) x(26)*ws_ts(x,i) + x(25);
vri = @(x,i) vei(x,i) - x(3);

tl = {'$v_{r1}$','$v_{r2}$', '$v_{r3}$'};

N = size(xk,2);
xk_vri = zeros(3,N);
xm_vri = zeros(3,N);

for k=1:N
    for i=0:2
        xk_vri(i+1,k) = vri(xk(:,k),i);
        xm_vri(i+1,k) = vri(xm(:,k),i);
    end
end

for i=1:3
    figure
    plot(t,xk_vri(i,:), t,xm_vri(i,:));
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$Velocity [\frac{m}{s}]$', 'Interpreter', 'latex', 'FontSize', 14);
    % grid on;
    legend('UKF', 'MPC');
    title(['$Estimation $' tl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end
end