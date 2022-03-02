%Define the MPC object with Yalmip:
 
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

Epsilon = Tau - Psi*P0 -Upsilon*Uprev;

P = Acal*P0 + Bcal_u*Uprev + Bcal_du*deltaU; % P: P*: Optimal states in prediction window

ident = eye(Lk);
for i=1:Hu-1
    ident = [ident; eye(Lk)];
end

% V = [ident ident];
% for i=1:Hu-1
%     temp = [zeros(i*Yk,Lk); ident(1:end-2*i,:)];
%     V = [V temp]; 
% end
% 
% deltaU_full = [Uprev; deltaU];
% U = V*deltaU_full; % U: u*: The optimal control window

%%%%%%%%%%%%
%%% Cost %%%
%%%%%%%%%%%%
Gcal = 2*Theta'*Qcal*Epsilon;
Hcal = Theta'*Qcal*Theta + Rcal;

Cost = Epsilon'*Qcal*Epsilon - deltaU'*Gcal + deltaU'*Hcal*deltaU;

%%%%%%%%%%%%%%%%%%%
%%% Constraints %%%
%%%%%%%%%%%%%%%%%%%

Constraints = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLACE YOUR CODE HERE (END)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Yalmip optimizer-object used for simulation and TurtleBot3 control
ops = sdpsettings('solver','sedumi');
MPCobj=optimizer(Constraints,Cost,ops,{P0,Uprev,refht},{U, P});
% U: u*: The optimal control window
% P: P*: Optimal states in prediction window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOU CAN PLACE YOUR OWN FUNCTIONS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


