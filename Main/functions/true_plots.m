function [] = true_plots(Lk,y,xk,xm,x_ul,x_vl,var_names,t)
tl = {'$P_e$','$v_r$'};

figure
plot(t,y(13,:)');
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Watt [W]$', 'Interpreter', 'latex', 'FontSize', 14);
% grid on;
legend('Non-linear');
title(['$Estimation $' tl(1)], 'Interpreter', 'latex', 'FontSize', 15);


figure
% plot(t,y(14,:), t,(xk(26,:)+xk(25,:)), t,(xm(26,:)+xm(25,:)));
plot(t,(xk(26,:)+xk(25,:)), t,(xm(26,:)+xm(25,:)));
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Velocity [\frac{m}{s}]$', 'Interpreter', 'latex', 'FontSize', 14);
% grid on;
legend('Non-linear', 'Linear');
title(['$Estimation $' tl(2)], 'Interpreter', 'latex', 'FontSize', 15);

filename = sprintf('plots/MPC_%s', var_names(1));
set(gcf, 'PaperPosition', [0 0 14.5 9.4]); %Position plot at left hand corner with width 25 and height 20.
set(gcf, 'PaperSize', [14.5 9.4]); %Set the paper to have width 25 and height 20.
saveas(gcf, filename, 'pdf') %Save figure

for i = 1:Lk
    figure
    % plot(t,xt(i,:),t,xk(i,:),t,xm(i,:));
    plot(t,xk(i,:),t,xm(i,:));
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    % grid on;
    legend('Non-linear', 'Linear');
    title(['$Estimation $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
    
    filename = sprintf('plots/MPC_%s', var_names(1+i));
    set(gcf, 'PaperPosition', [0 0 14.5 9.4]); %Position plot at left hand corner with width 25 and height 20.
    set(gcf, 'PaperSize', [14.5 9.4]); %Set the paper to have width 25 and height 20.
    saveas(gcf, filename, 'pdf') %Save figure
end
end

