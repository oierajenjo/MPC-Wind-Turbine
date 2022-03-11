%Define the MPC object with Yalmip:
clear all
clc
variables_IPC
% Constraints limit values
U_c.theta_min = 0;
U_c.theta_max = pi/2;
U_c.Tg_min = 1;
U_c.Tg_max = 2.159e7;

dU_c.theta_dot = deg2rad(9); % Pitch angle max and min angular speed

Z_c.omega_r_min = convangvel(5,'rpm', 'rad/s');
Z_c.omega_r_max = - convangvel(10,'rpm', 'rad/s');
Z_c.xtdd_min = 0; 
Z_c.xtdd_max = -0;
Z_c.ytdd_min = 0;
Z_c.ytdd_max = -0;
Z_c.My_min1 = B.m*B.l*B.xdd_min; Z_c.My_max1 = -B.m*B.l*B.xdd_max;
Z_c.My_min2 = B.m*B.l*B.xdd_min; Z_c.My_max2 = -B.m*B.l*B.xdd_max;
Z_c.My_min3 = B.m*B.l*B.xdd_min; Z_c.My_max3 = -B.m*B.l*B.xdd_max;
Z_c.Mx_min1 = B.m*B.l*B.ydd_min; Z_c.Mx_max1 = -B.m*B.l*B.ydd_max;
Z_c.Mx_min2 = B.m*B.l*B.ydd_min; Z_c.Mx_max2 = -B.m*B.l*B.ydd_max;
Z_c.Mx_min3 = B.m*B.l*B.ydd_min; Z_c.Mx_max3 = -B.m*B.l*B.ydd_max;
Z_c.Pe_min = D.eta*U_c.Tg_min*Z_c.omega_r_min;
Z_c.Pe_max = -D.eta*U_c.Tg_max*Z_c.omega_r_max;
Z_c.vr_min = 4;
Z_c.vr_max = -25;
Z_c.azim_min = 0;
Z_c.azim_max = -0;

% Ts = Ts; % Sample time for the MPC controller
Ts = 0.05;
Hp = 10; % Prediction Horizon
Hu = 10; % Control Horizon
Hw = 1; % Window parameter

Lk = 27;
Uk = 4;
Yk = 12;

refht=sdpvar(Yk,Hp);      % P_{ref}: The window containing the pos-reference
P0=sdpvar(Lk,1);          % P(k):    The current state
Uprev=sdpvar(Uk,1);       % U(k-1):  The previous input command.
deltaU = sdpvar(Uk*Hu,1); % DeltaU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLACE YOUR CODE HERE (START)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% System Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%
Ampc = eye(Lk); % State Matrix

Bmpc = Ts*eye(Lk,Uk); % Input Matrix

Cmpc = eye(Yk,Lk);

Q_c = 1; % Error Weight
Qmpc = Q_c*eye(Yk);

R_c = 0; % Input Weight
Rmpc = R_c*eye(Uk);

% R_d = 10; % Change Input Weight
% R_del = R_d*eye(Uk);

%%%%%%%%%%%%%%%
%%% Lifting %%%
%%%%%%%%%%%%%%%
Acal = [];
Ai = Ampc;
Bcal_u = [];
Bi = zeros(Lk,Uk);
for i=1:Hp
    if i~=1
        Ai = Ai*Ampc;
    end
    Acal = [Acal; Ai];
    
    Bi = Ampc*Bi+Bmpc;
    Bcal_u = [Bcal_u; Bi];
end

Bcal_du = Bcal_u;
for i=1:Hu-1
    temp = [zeros(i*Lk,Uk); Bcal_u(1:end-Lk*i,:)];
    Bcal_du = [Bcal_du temp];
end

Ccal = repmat({Cmpc}, 1, Hp);
Ccal = blkdiag(Ccal{:});
Ccal = Ccal(Yk*(Hw-1)+1:end,:);

Psi = Ccal*Acal;
Upsilon = Ccal*Bcal_u;
Theta = Ccal*Bcal_du;

Tau = reshape(refht,[Yk*Hp 1]); % Convert to column vector
Tau = Tau(Yk*(Hw-1)+1:end,:);

Qcal = repmat({Qmpc}, 1, Hp);
Qcal = blkdiag(Qcal{:});
Qcal = Qcal(Yk*(Hw-1)+1:end,:);

Rcal = repmat({Rmpc}, 1, Hu);
Rcal = blkdiag(Rcal{:});
% Rcal_del = repmat({R_del}, 1, Hu);
% Rcal_del = blkdiag(Rcal_del{:});

Epsilon = Tau - Psi*P0 - Upsilon*Uprev;

P = Acal*P0 + Bcal_u*Uprev + Bcal_du*deltaU; % P: P*: Optimal states in prediction window

ident = eye(Uk);
for i=1:Hu-1
    ident = [ident; eye(Uk)];
end

deltaU_full = [Uprev; deltaU];
V = [ident ident];
for i=1:Hu-1
    temp = [zeros(i*Uk,Uk); ident(1:end-Uk*i,:)];
    V = [V temp]; 
end

deltaU_full = [Uprev; deltaU];
U = V*deltaU_full; % U: u*: The optimal control window

%%%%%%%%%%%%
%%% Cost %%%
%%%%%%%%%%%%
Gcal = 2*Theta'*Qcal*Epsilon;
Hcal = Theta'*Qcal*Theta + Rcal;

Cost = Epsilon'*Qcal*Epsilon - deltaU'*Gcal + deltaU'*Hcal*deltaU;

%%%%%%%%%%%%%%%%%%%
%%% Constraints %%%
%%%%%%%%%%%%%%%%%%%
% Actuator Range Constraints
f = [-1; 1];
f_rep = repmat({f}, 1, Uk*Hu);
F_L = blkdiag(f_rep{:});
u_cons = [U_c.theta_min; -U_c.theta_max; U_c.theta_min; -U_c.theta_max; 
        U_c.theta_min; -U_c.theta_max; U_c.Tg_min; -U_c.Tg_max];
U_L = u_cons;
for i=1:Hu-1
    U_L = [U_L; u_cons]; 
end
F = [F_L U_L];

% Actuator Slew Rate Constraints
e = [-1; 1];
e_rep = repmat({e}, 1, Uk-1);
e_rep{end+1} = zeros(2,1);
e_L = blkdiag(e_rep{:});

e_L = repmat({e_L}, 1, Hu);
E_L = blkdiag(e_L{:});
du_cons = [-dU_c.theta_dot*ones(6,1); zeros(2,1)];
dU_L = du_cons;
for i=1:Hu-1
   dU_L = [dU_L; du_cons]; 
end
E = [E_L dU_L];

% Actuator Range Constraints
g = [-1; 1];
g_rep = repmat({g}, 1, Yk-1);
g_rep{end+1} = [-1;0];
g_L = blkdiag(g_rep{:});

g_L = repmat({g_L}, 1, Hu);
G_L = blkdiag(g_L{:});
z_cons = cell2mat(struct2cell(Z_c));
Z_L = z_cons;
for i=1:Hu-1
    Z_L = [Z_L; z_cons]; 
end
G = [G_L Z_L];
Zcal = Psi*P0 + Upsilon*Uprev + Theta*deltaU;

Constraints = [F*[U;1]<=0; E*[deltaU;1]<=0; G*[Zcal;1]<=0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLACE YOUR CODE HERE (END)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Yalmip optimizer-object used for simulation and TurtleBot3 control
ops = sdpsettings('solver','sedumi');
MPCobj=optimizer(Constraints,Cost,ops,{P0,Uprev,refht},{U,P});
% U: u*: The optimal control window
% P: P*: Optimal states in prediction window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOU CAN PLACE YOUR OWN FUNCTIONS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


