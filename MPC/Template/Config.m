%Initialization file for some Parameters: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION (You can play around with the configuration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulation/Implementation time:
T = 30; % Final time
Tsim = 0.01; %Simulation time step:

%Initial state
X0 = [0;0;0];

% Parameters for the reference trajectory
omega = 0.25;
radius = 0.5;
Phi = pi/4;
x_offset = -0.8;
y_offset = -0.4;

%Define the limit for the state constraint
y_constraint=-0.6;