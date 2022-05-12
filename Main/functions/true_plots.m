function [] = true_plots(Lk,y,xk,xt,xm,x_ul,x_vl,var_names,t)
     tl = {'$P_e$','$v_r$'};
%     figure
%     plot(t,y(13,:)');
%     xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
%     ylabel('$Watt [W]$', 'Interpreter', 'latex', 'FontSize', 14);
%     grid on;
%     legend('True');
%     title(['$Estimation $' tl(1)], 'Interpreter', 'latex', 'FontSize', 15); 
    
    filename = sprintf('MPC_%s', var_names(1));

    figure(1)
    %plot(t,y(14,:), t,(xk(26,:)+xk(25,:)), t,(xm(26,:)+xm(25,:)));
    plot(t,(xk(26,:)+xk(25,:)), t,(xm(26,:)+xm(25,:)));
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$Velocity [\frac{m}{s}]$', 'Interpreter', 'latex', 'FontSize', 14);
%     grid on;
    legend('UKF', 'MPC');
    title(['$Estimation $' tl(2)], 'Interpreter', 'latex', 'FontSize', 15);
    set(gcf, 'PaperPosition', [0 0 6 3.7]); %Position plot at left hand corner with width 25 and height 20.
    set(gcf, 'PaperSize', [5.7 3.7]); %Set the paper to have width 25 and height 20.
    saveas(gcf, filename, 'pdf') %Save figure
    % [0 0 13 8]
    % [12 8]
for i = 1:Lk
    figure(i+1)
    %plot(t,xt(i,:),t,xk(i,:),t,xm(i,:));
    plot(t,xk(i,:),t,xm(i,:));
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel(x_ul(i), 'Interpreter', 'latex', 'FontSize', 14);
%     grid on;
    legend('UKF', 'MPC');
    title(['$Estimation $' x_vl(i)], 'Interpreter', 'latex', 'FontSize', 15);
    
    filename = sprintf('MPC_%s', var_names(1+i));

    set(gcf, 'PaperPosition', [0 0 6 3.7]); %Position plot at left hand corner with width 25 and height 20.
    set(gcf, 'PaperSize', [5.7 3.7]); %Set the paper to have width 25 and height 20.
    saveas(gcf, filename, 'pdf') %Save figure
end
end

