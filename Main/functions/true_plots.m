function [] = true_plots(Lk,yt,xt,x_ul,x_vl,t)
    pe = {'$P_e$'};
    figure
    plot(t,yt(10,:)');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$Watt [W]$', 'Interpreter', 'latex', 'FontSize', 14);
%     grid on;
    legend('True');
    title(['$Estimation $' pe(1)], 'Interpreter', 'latex', 'FontSize', 15);

for i = 1:Lk
    figure
    plot(t,xt(i,:));
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
%     grid on;
    legend('UKF');
    title(['$Estimation $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end
end

