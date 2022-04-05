function [] = true_plots(yt,y_me,xt,data,t)
    figure
    plot(t,yt(1,:)',t,y_me(1,:));
    legend('yt','Bladed')
    % ylim([-2.6 2.6]);
    title('wr');
    xlabel('s')

    figure
    plot(t,yt(2,:)',t,y_me(2,:));
    legend('yt','Bladed')
    % ylim([-2.6 2.6]);
    title('xtddot');
    xlabel('s')

    figure
    plot(t,yt(3,:)',t,y_me(3,:));
    legend('yt','Bladed')
    % ylim([-2.6 2.6]);
    title('ytddot');
    xlabel('s')

    figure
    plot(t,yt(4,:)',t,y_me(4,:));
    legend('yt','Bladed')
    % ylim([-2.6 2.6]);
    title('My1');
    xlabel('s')

    figure
    plot(t,yt(7,:)',t,y_me(7,:));
    legend('yt','Bladed')
    % ylim([-2.6 2.6]);
    title('Mx1');
    xlabel('s')

    figure
    plot(t,yt(10,:)',t,y_me(10,:));
    legend('yt','Bladed')
    % ylim([-2.6 2.6]);
    title('Pe');
    xlabel('s')

    figure
    plot(t,yt(11,:)',t,y_me(11,:));
    legend('yt','Bladed')
    % ylim([-2.6 2.6]);
    title('vr');
    xlabel('s')

    figure
    plot(t,xt(2,:)',t,data.Data(:,224));
    legend('xt','Bladed')
    % ylim([-2.6 2.6]);
    title('xt');
    xlabel('s')

    figure
    plot(t,xt(4,:)',t,data.Data(:,225));
    legend('xt','Bladed')
    % ylim([-2.6 2.6]);
    title('yt');
    xlabel('s')

    figure
    plot(t,xt(6,:)',t,data.Data(:,85));
    legend('xt','Bladed')
    % ylim([-2.6 2.6]);
    title('xb1');
    xlabel('s')

    figure
    plot(t,xt(12,:)',t,data.Data(:,86));
    legend('xt','Bladed')
    % ylim([-2.6 2.6]);
    title('yb1');
    xlabel('s')

    figure
    plot(t,xt(24,:)',t,data.Data(:,20));
    legend('xt','Bladed')
    % ylim([-2.6 2.6]);
    title('Tg');
    xlabel('s')
end

