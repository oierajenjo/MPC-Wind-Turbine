%Initialization file for some Parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONFIGURATION (You can play around with the configuration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
variables_IPC
%Simulation/Implementation time:
T = N*Ts; % Final time
Tsim = Ts; %Simulation time step:

%Initial state
X0 = x_i;

% Parameters for the reference trajectory

%Define the limit for the state constraint
