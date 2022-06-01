close all
clc

yt(:,1) = y_me(:,1);
Pe_f = @(x) D.eta*x(24)*x(1);
for k=1:N
    pe(k) = Pe_f(x_tv(:,k));
    lam(k) = lamb_eq(x_tv(:,k));
end
tl = {'$P_e$','$\lambda$'};

figure
plot(t,z_mpc(end,:),t,pe);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$Watt [W]$', 'Interpreter', 'latex', 'FontSize', 14);
% grid on;
legend('Non-linear','Expected');
title(['$Estimation $' tl(1)], 'Interpreter', 'latex', 'FontSize', 15);
filename = 'plots/MPC_Pe';
set(gcf, 'PaperPosition', [0 0 14.5 9.4]); %Position plot at left hand corner with width 25 and height 20.
set(gcf, 'PaperSize', [14.5 9.4]); %Set the paper to have width 25 and height 20.
saveas(gcf, filename, 'pdf') %Save figure

figure
plot(t,z_mpc(10,:)',t,lam);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\cdot$', 'Interpreter', 'latex', 'FontSize', 14);
% grid on;
legend('Non-linear','Expected');
title(['$Estimation $' tl(2)], 'Interpreter', 'latex', 'FontSize', 15);
filename = 'plots/MPC_lambda';
set(gcf, 'PaperPosition', [0 0 14.5 9.4]); %Position plot at left hand corner with width 25 and height 20.
set(gcf, 'PaperSize', [14.5 9.4]); %Set the paper to have width 25 and height 20.
saveas(gcf, filename, 'pdf') %Save figure