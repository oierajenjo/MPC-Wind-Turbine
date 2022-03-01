% ROS interface with Matlab
% initilize script with some parameters: 
%exampleHelperROSShutDownSampleNetwork
rosshutdown
clear
%close all 
%Init
Config
MPCdefinition;
%%%% 
% Reference
%%%%
ref = [x_offset+radius*cos(omega*(0:Ts:(T+Hp))+Phi);...
       y_offset+radius*sin(omega*(0:Ts:(T+Hp))+Phi)]; % Position

dref = [ omega*radius*sin(omega*(0:Ts:(T+Hp))+Phi);...
        -omega*radius*cos(omega*(0:Ts:(T+Hp))+Phi)];
% Connect to the ROS on the robot:
rosinit('192.168.1.101')
% start a publisher and a message template on the reset topic 
[pub,msg] = rospublisher('reset','std_msgs/Empty')
% send a message to reset odometry 
send(pub,msg)
%Subscribe to odometry:
states= rossubscriber('/odom');
%Creat a publisher for heading with a heading msg:
headingpub = rospublisher('/heading', 'geometry_msgs/Vector3');
headingmsg = rosmessage(headingpub);
%Create a publisher for velocity with a velocity msg:
velpub = rospublisher('/cmd_vel_x', 'geometry_msgs/Vector3');
velmsg = rosmessage(velpub);
%Start the MPC:  %Make sure your robot is on the zero position. if you turn
%on and don't move it from its position then it will be there.
time=0:Ts:T;
%initialize vectors:
X=zeros(3,length(time)); %States [x;y;th] th is the heading angle
V=zeros(1,length(time)); %Velocity magnitude (Velocity in Body frame)
W=zeros(1,length(time)); %Angular Velocity
vx=zeros(1,length(time)); %MPC output of velocity on x
vy=zeros(1,length(time)); %MPC output of velocity on y
thd=zeros(1,length(time)); %Desired heading angle from MPC
%[ysound, Fs] =audioread('nn.mp3');
%sound(ysound, Fs, 16);
uprev=[0;0];
pause(3)
headingmsg.Z=0;
send(headingpub,headingmsg)

for k=1:length(time)
Data=receive(states,1);
%Current x,y, and th
X(1,k)=Data.Pose.Pose.Position.X;
X(2,k)=Data.Pose.Pose.Position.Y;
%To read the heading: 
qx=Data.Pose.Pose.Orientation.X;
qy=Data.Pose.Pose.Orientation.Y;
qz=Data.Pose.Pose.Orientation.Z;
qw=Data.Pose.Pose.Orientation.W;
Euler=quat2eul([qx,qy,qz,qw]);
X(3,k)=wrapToPi(Euler(3));
%Apply the MPC:
tic
res=MPCobj({X(1:2,k),uprev,ref(:,k:k+(Hp-1))});
Comptime=toc;
if(Comptime>0.5)
    disp('large time')
    disp(Comptime)
end
u_L = res{1};
u = reshape(u_L, [2, Hu]);
% u_pre(:,:,k) = u;
% P_L = res{2};
% P_pre(:,:,k) = reshape(P_L, [2, Hp]);
uprev=u(:,1);
%u=MPCcode(ref(:,k:k+(H-1)),X(1:2,k),Ts,H); %Pass a reference that matches the sample time of
%MPC.
vx(k)=u(1,1);
vy(k)=u(2,1);
V(k)=sqrt(vx(k)^2+vy(k)^2);
%Desired heading angle: 
thd(k)=atan2(vy(k),vx(k));
headingmsg.Z=thd(k);
velmsg.X=V(k);
send(headingpub,headingmsg)
send(velpub,velmsg);
pause(Ts-Comptime)
end
%To insure the robot is stopping
velmsg.X=0;
send(velpub,velmsg);
%exampleHelperROSShutDownSampleNetwork

%%
figure
rosshutdown
h1=plot(X(1,:),X(2,:));
hold on 
%h4=plot(X(1,1:Tsim/Ts:end),X(2,1:Tsim/Ts:end),'bX');
%hold on 
h2=plot(ref(1,:),ref(2,:),'g');
hold on 
h3=plot(ref(1,:),-0.6*ones(1,length(ref(1,:))),'r');
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
title('MPC implementation','interpreter','latex')
grid on 
leg=legend('Actualt Trajectory',...
'Desired Trajectory','Constraint');
set(leg,'interpreter','latex');
