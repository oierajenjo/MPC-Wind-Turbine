% MPCsim

% Description: This file simulates the non-linear model

clear all; close all; clc

addpath('functions/')

%Run the config script for reference generation:
Config

%Define the MPC controller:
MPCdefinition;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate reference trajectory: (Circle)
ref = [x_offset+radius*cos(omega*(0:Ts:(T+Hp))+Phi);...
       y_offset+radius*sin(omega*(0:Ts:(T+Hp))+Phi)]; % Position

dref = [ omega*radius*sin(omega*(0:Ts:(T+Hp))+Phi);...
        -omega*radius*cos(omega*(0:Ts:(T+Hp))+Phi)]; % Velocity

%Time vector for simulation:
tsim=0:Tsim:T;

%initialize vectors for simulation:
X = zeros(3,length(tsim)); %States [x;y;th] th is the heading angle
V = zeros(1,length(tsim)); %Velocity magnitude (Velocity in Body frame)
W = zeros(1,length(tsim)); %Angular Velocity
vx = zeros(1,length(tsim)); %MPC output of velocity on x
vy = zeros(1,length(tsim)); %MPC output of velocity on y
thd = zeros(1,length(tsim)); %Desired heading angle from MPC

P_pre = cell(length(tsim),1); %Container for prediction window of states
u_pre = cell(length(tsim),1); %Container for prediction window of inputs

%Initial states: 
X(:,1)=X0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indref=1; %index for reference. increases everytime the MPC is excuted
uprev=[0;0];
for k=1:length(tsim)-1
    if(mod(tsim(k),Ts)==0) %Is it time to excute the MPC? 
        res=MPCobj({X(1:2,k),uprev,ref(:,indref+1:indref+Hp)});
        
        u_L = res{1};
        u = reshape(u_L, [2, length(u_L)/2]);
        u_pre{k} = u;
        uprev=u(:,1);
        
        P_L = res{2};
        P_pre{k} = reshape(P_L, [2, length(P_L)/2]);
        
        %MPC.
        vx(k)=u(1,1);
        vy(k)=u(2,1);
        V(k)=sqrt(vx(k)^2+vy(k)^2);
        
        %Desired heading angle: 
        thd(k)=atan2(vy(k),vx(k));
        
        % Update index
        indref=indref+1;
    else
        vx(k)=vx(k-1);
        vy(k)=vy(k-1);
        V(k)=V(k-1);
        thd(k)=thd(k-1);
    end
    
    %P controller for heading:
    diff = wrapToPi(thd(k)-X(3,k));
    W(k)=2.8*(diff);
    
    %Apply model: 
    X(1,k+1)=X(1,k)+V(k)*cos(X(3,k))*Tsim;
    X(2,k+1)=X(2,k)+V(k)*sin(X(3,k))*Tsim;
    X(3,k+1)=X(3,k)+W(k)*Tsim;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure

corder = get(gca,'colororder');

%Post plotting
h1=plot(X(1,:),X(2,:));
hold on 
% Just calculating and setting figure limits
idx = 1;
buffer = 0.1*(max(max(ref(idx,:)), max(X(idx,:)))-min(min(ref(idx,:)), min(X(idx,:))));
xlimit = [min(min(ref(idx,:)), min(X(idx,:)))-buffer max(max(ref(idx,:)), max(X(idx,:)))+buffer];
set(gca, 'Xlim', xlimit)

idx = 2;
buffer = 0.1*(max(max(ref(idx,:)), max(X(idx,:)))-min(min(ref(idx,:)), min(X(idx,:))));
ylimit = [min(min(ref(idx,:)), min(X(idx,:)))-buffer max(max(ref(idx,:)), max(X(idx,:)))+2*buffer];
set(gca, 'Ylim', ylimit)

h4=plot(X(1,1:round(Ts/Tsim):end),X(2,1:round(Ts/Tsim):end),'bX');
h2=plot(ref(1,:),ref(2,:),'g');
h3=line(xlimit,[y_constraint y_constraint], 'color', 'r');
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
title('MPC simulation','interpreter','latex')
grid on 
leg=legend('Actualt Trajectory','Trajectory at MPC Sample Time',...
'Desired Trajectory','Constraint');
set(leg,'interpreter','latex');


% Velocity and heading plot
figure 
% Vx
subplot(3,1,1)
plot(tsim,vx)
xlabel('Time [s]','interpreter','latex')
ylabel('Velocity $v_{x}$ [m/s]','interpreter','latex')
leg=legend('$v_{x}$');
set(leg,'interpreter','latex');
title('Velocity along x-axis (Prediction Model)')
grid on 

% Vy
subplot(3,1,2)
plot(tsim,vy)
xlabel('Time [s]','interpreter','latex')
ylabel('Velocity $v_{x}$ [m/s]','interpreter','latex')
leg=legend('$v_{y}$');
set(leg,'interpreter','latex');
title('Velocity along y-axis (Prediction Model)')
grid on 

% Heading
subplot(3,1,3)
plot(tsim,wrapToPi(X(3,:)))
xlabel('Time [s]')
ylabel('Heading [rad]','interpreter','latex')
leg=legend('$\Theta_v$');
set(leg,'interpreter','latex');
title('Heading (Prediction Model)')
grid on 

% Animation of trajectory

figure('Position', [100, 100, 600, 600])
sp1 = subplot(2,1,1);

hold on


% Just calculating and setting figure limits
idx = 1;
buffer = 0.1*(max(max(ref(idx,:)), max(X(idx,:)))-min(min(ref(idx,:)), min(X(idx,:))));
xlimit = [min(min(ref(idx,:)), min(X(idx,:)))-buffer max(max(ref(idx,:)), max(X(idx,:)))+buffer];
set(gca, 'Xlim', xlimit)

idx = 2;
buffer = 0.1*(max(max(ref(idx,:)), max(X(idx,:)))-min(min(ref(idx,:)), min(X(idx,:))));
ylimit = [min(min(ref(idx,:)), min(X(idx,:)))-buffer max(max(ref(idx,:)), max(X(idx,:)))+buffer];
set(gca, 'Ylim', ylimit)
grid on

% References
href=plot(ref(1,:),ref(2,:),'g');

% y constraint
line(xlimit,[y_constraint y_constraint], 'color', 'r')

% cart
L = 0.15;
W = 0.1;
cart = ani_cart(L, W , [0.4940, 0.1840, 0.5560], 1.2);
cart.plot_cart(X(1:2,1), V(1), wrapToPi(X(3,1)));

P_pre_1 = P{1};

% Prediction horizon
ref_x = ref(1,1:Hp);
ref_y = ref(2,1:Hp);
ref_curve = ani_curve(corder(5,:),2);
ref_curve.plot_curve(ref_x, ref_y)

% Prediction horizon
P_pre_x = P_pre_1(1,:);
P_pre_y = P_pre_1(2,:);
curve = ani_curve([0.8500, 0.3250, 0.0980],1.5);
curve.plot_curve(P_pre_x, P_pre_y)

% Path
x_path = X(1,1);
y_path = X(2,1);
path = ani_curve([0.9290, 0.6940, 0.1250], 1);
path.plot_curve(x_path, y_path)


title('Trajectory Animation', 'fontsize', 14)
xlabel('Pos. x', 'fontsize', 12)
ylabel('Pos. y', 'fontsize', 12)
daspect([1 1 1])
pos1 = get(sp1, 'Position'); % gives the position of current sub-plot
resize = 0.10;
new_pos1 = pos1+[-resize -1.2*resize 1.5*resize 1.5*resize];
set(sp1, 'Position',new_pos1);
sp1_pos = get(sp1,  'Position');


sp2 = subplot('Position', [new_pos1(1)+0.08,new_pos1(2)-0.38,new_pos1(3)-0.155,new_pos1(4)-0.25]);
hold on
xlabel('Time [s]', 'fontsize', 12)
ylabel('u [m/s]', 'fontsize', 12)
set(gca, 'xlim', [0, T])
buffer = 0.1;
vlim = 0.22;
ylim = vlim/sqrt(2);
set(gca, 'ylim', [-ylim-buffer, ylim+buffer])
title('Input', 'fontsize', 14)

vx_curve = vx(1);
vy_curve = vy(1);
v_curve = sqrt(vx(1).^2+vy(1).^2);
t_curve = tsim(1);
plt1 = plot(t_curve, vx_curve, 'linewidth', 1.5);
plt2 = plot(t_curve, vy_curve, 'linewidth', 1.5);
plt3 = plot(t_curve, v_curve, 'linewidth', 1.5);

line([0, T],[-ylim -ylim],'color', 'k', 'linestyle','--', 'linewidth', 1)
line([0, T],[ylim ylim],'color', 'k', 'linestyle','--', 'linewidth', 1)
line([0, T],[-vlim -vlim],'color', [0.9290    0.6940    0.1250], 'linestyle','--', 'linewidth', 1)
line([0, T],[vlim vlim],'color', [0.9290    0.6940    0.1250], 'linestyle','--', 'linewidth', 1)

plt1.XDataSource = 't_curve';
plt1.YDataSource = 'vx_curve';

plt2.XDataSource = 't_curve';
plt2.YDataSource = 'vy_curve';

plt3.XDataSource = 't_curve';
plt3.YDataSource = 'v_curve';

leg = legend('$v_x$', '$v_y$', '$v$');
set(leg, 'location', 'northeast');
set(leg,'interpreter','latex');
grid on
indref=1;
for k = 1:length(tsim)-1
    P_pre_k = P_pre{k};
    if(mod(tsim(k),Ts)==0) %Is it time to excute the MPC? 
        ref_curve.update_curve(ref(1,indref+1:indref+Hp), ref(2,indref+1:indref+Hp))
        curve.update_curve(P_pre_k(1,:), P_pre_k(2,:))
        indref = indref + 1;
    end
    
    t_curve = tsim(1:k);
    vx_curve = vx(1:k);
    vy_curve = vy(1:k);
    v_curve = sqrt(vx(1:k).^2+vy(1:k).^2);
    
    set(plt1, 'xdata', t_curve);
    set(plt1, 'ydata', vx_curve);
    set(plt2, 'xdata', t_curve);
    set(plt2, 'ydata', vy_curve);
    set(plt3, 'xdata', t_curve);
    set(plt3, 'ydata', v_curve);
    
    cart.update_cart(X(1:2,k), V(k), wrapToPi(X(3,k)));
    path.update_curve(X(1,1:k), X(2,1:k));
    pause(0.001)
end


