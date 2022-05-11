function [] = true_plots(Lk,y,xk,xt,xm,x_ul,x_vl,t)
tl = {'$P_e$','$v_r$'};
figure
plot(t,y(13,:)');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Watt [W]$', 'Interpreter', 'latex', 'FontSize', 14);
%     grid on;
legend('True');
title(['$Estimation $' tl(1)], 'Interpreter', 'latex', 'FontSize', 15);

figure
plot(t,y(14,:), t,(xk(26,:)+xk(25,:)), t,(xm(26,:)+xm(25,:)));
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Velocity [\frac{m}{s}]$', 'Interpreter', 'latex', 'FontSize', 14);
%     grid on;
legend('True', 'UKF', 'MPC');
title(['$Estimation $' tl(2)], 'Interpreter', 'latex', 'FontSize', 15);

for i = 1:Lk
    figure
    plot(t,xt(i,:),t,xk(i,:),t,xm(i,:));
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    %     grid on;
    legend('True', 'UKF', 'MPC');
    title(['$Estimation $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end
end

