%Define the MPC object with Yalmip:

Ts = 0.5; % Sample time for the MPC controller

Hp = 10; % Prediction Horizon

Hu = 10; % Control Horizon

refht = sdpvar(2,Hp);     % P_{ref}: The window containing the pos-reference
P0 = sdpvar(2,1);         % P(k):    The current state
Uprev = sdpvar(2,1);      % U(k-1):  The previous input command.
deltaU = sdpvar(2*Hu,1);  % DeltaU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLACE YOUR CODE HERE (START)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% System Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%
A = eye(2); % State Matrix
n = size(A,1);

B = Ts*eye(2); % Input Matrix
m = size(B,1);

Q_c = 1; % Error Weight
Q = Q_c*eye(n);

R_c = 0; % Input Weight
R = R_c*eye(m); 

R_d = 10; % Change Input Weight
R_del = R_d*eye(m);

%%%%%%%%%%%%%%%
%%% Lifting %%%
%%%%%%%%%%%%%%%
Acal = [];
Ai = A;
Bcal_u = [];
Bi = zeros(n,m);
for i=1:Hp
    if i~=1
        Ai = Ai*A;
    end
    Acal = [Acal; Ai];
    
    Bi = A*Bi+B;
    Bcal_u = [Bcal_u; Bi];
end

Bcal_du = Bcal_u;
for i=1:Hu-1
    temp = [zeros(i*m,n); Bcal_u(1:end-2*i,:)];
    Bcal_du = [Bcal_du temp];
end

Qcal = repmat({Q}, 1, Hp);
Qcal = blkdiag(Qcal{:});
Rcal = repmat({R}, 1, Hu);
Rcal = blkdiag(Rcal{:});
Rcal_del = repmat({R_del}, 1, Hu);
Rcal_del = blkdiag(Rcal_del{:});

P = Acal*P0 + Bcal_u*Uprev + Bcal_du*deltaU; % P: P*: Optimal states in prediction window

ident = eye(n);
for i=1:Hu-1
    ident = [ident; eye(n)];
end

V = [ident ident];
for i=1:Hu-1
    temp = [zeros(i*m,n); ident(1:end-2*i,:)];
    V = [V temp]; 
end

deltaU_full = [Uprev; deltaU];
U = V*deltaU_full; % U: u*: The optimal control window

%%%%%%%%%%%%
%%% Cost %%%
%%%%%%%%%%%%

refht_col = reshape(refht,[2*size(refht,2) 1]); % Convert to column vector

Cost = Ts*((P-refht_col)'*Qcal*(P-refht_col) + deltaU'*Rcal_del*deltaU +...
       U'*Rcal*U);

%%%%%%%%%%%%%%%%%%%
%%% Constraints %%%
%%%%%%%%%%%%%%%%%%%
g = [0 -1];
g_rep = repmat({g}, [1 Hp]);
G_L = blkdiag(g_rep{:});
GL = [G_L ones(Hp,1)*y_constraint*0.95];%(y_constraint+0.15/2)*0.99];
    
max_v = 0.22/sqrt(n);

Constraints = [GL*[P;1]<=0,-max_v<=U, U<=max_v];

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


