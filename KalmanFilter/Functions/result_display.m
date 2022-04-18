function result_display(t, Lk, xk, xt, x_me, x_ul, x_vl)
for i = 1:Lk
    figure
    subplot(1,2,1);
    %     plot(t,xt(i,:),'r-','LineWidth', 2);
    plot(t,xk(i,:),'b-', t,xt(i,:),'r-', t,x_me(i,:),'g-');
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF','True','Bladed');
    title(['$Estimations $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
    
    subplot(1,2,2);
    plot(t,xk(i,:)-xt(i,:),'b-');
    xlabel('$Time (s)$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    legend('UKF');
    title(['$Error $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
end
end

