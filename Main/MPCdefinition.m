%Define the MPC object with Yalmip:
refht=sdpvar(Zk,Hp);      % Z_{ref}: The window containing the pos-reference
X0=sdpvar(Lk,1);          % X(k):    The current state
Uprev=sdpvar(Uk,1);       % U(k-1):  The previous input command.
deltaU = sdpvar(Uk*Hu,1); % DeltaU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLACE YOUR CODE HERE (START)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%% System Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%

A1 = [zeros(1,11), -a1*ones(1,3);
    zeros(1,2), 1, zeros(1,11);
    0, -b3, -b4, zeros(1,2), b1*ones(1,3), b2*ones(1,3), zeros(1,3);
    zeros(1,4), 1, zeros(1,9);
    zeros(1,3), -c4, -c5, zeros(1,6), c2*ones(1,3);
    zeros(3,8), eye(3), zeros(3,3);
    d10(xeq,0), d3, d4(xeq,0), zeros(1,2), -d1, zeros(1,2), d2(xeq,0), zeros(1,5);
    d10(xeq,1), d3, d4(xeq,1), zeros(1,3), -d1, zeros(1,2), d2(xeq,1), zeros(1,4);
    d10(xeq,2), d3, d4(xeq,2), zeros(1,4), -d1, zeros(1,2), d2(xeq,2), zeros(1,3);
    zeros(3,14)];

A2 = [zeros(1,9), -a2, zeros(1,3);
    zeros(3,13);
    c3*ones(1,3), zeros(1,6), -c1, zeros(1,3);
    zeros(3,13);
    d8(xeq,0), zeros(1,2), d7(xeq,0), zeros(1,6), d6(xeq,0), d5(xeq,0), d9(xeq,0);
    0, d8(xeq,1), zeros(1,2), d7(xeq,1), zeros(1,5), d6(xeq,1), d5(xeq,1), d9(xeq,1);
    0, 0, d8(xeq,2), zeros(1,2), d7(xeq,2), zeros(1,4), d6(xeq,2), d5(xeq,2), d9(xeq,2);
    eye(3), zeros(3,10)];

A3 = [-e9(xeq,0), 0, -e7(xeq,0), e3, e4, zeros(1,3), -e8(xeq,0), 0, 0, -e1, 0, 0;
    -e9(xeq,1), 0, -e7(xeq,1), e3, e4, zeros(1,4), -e8(xeq,1), 0, 0, -e1, 0;
    -e9(xeq,2), 0, -e7(xeq,2), e3, e4, zeros(1,5), -e8(xeq,2), 0, 0, -e1;
    zeros(9,14);
    1, zeros(1,13)];

A4 = [-e2(xeq,0), zeros(1,2), -e11(xeq,0), zeros(1,6), -e6(xeq,0), -e5(xeq,0), -e10(xeq,0);
    0, -e2(xeq,1), zeros(1,2), -e11(xeq,1), zeros(1,5), -e6(xeq,1), -e5(xeq,1), -e10(xeq,1);
    0, 0, -e2(xeq,2), zeros(1,2), -e11(xeq,2), zeros(1,4), -e6(xeq,2), -e5(xeq,2), -e10(xeq,2);
    zeros(3,6), eye(3), zeros(3,4);
    zeros(3), -f1*eye(3), -f2*eye(3), zeros(3,4);
    zeros(1,9), -g1, zeros(1,3);
    zeros(1,10), -W.w_p, zeros(1,2);
    zeros(2,13)];

Ampc = eye(Lk)+ Ts*[A1 A2; A3 A4];
Ampc = Sx\Ampc*Sx;% State Matrix


Bmpc = Ts*[zeros(3,Lk-7), f1*eye(3), zeros(3,4);
    zeros(1,Lk-4), g1, zeros(1,3)]';
Bmpc = Sx\Bmpc; % Input Matrix


lambda_rows = [
    r1(xeq,1), 0, r6(xeq,1), zeros(1,5), r7(xeq,1), zeros(1,5), r2(xeq,1), zeros(1,9), r5(xeq,1), r3(xeq,1), r4(xeq,1);
    r1(xeq,2), 0, r6(xeq,2), zeros(1,6), r7(xeq,2), zeros(1,5), r2(xeq,2), zeros(1,8), r5(xeq,2), r3(xeq,1), r4(xeq,1);
    r1(xeq,3), 0, r6(xeq,3), zeros(1,7), r7(xeq,3), zeros(1,5), r2(xeq,3), zeros(1,7), r5(xeq,3), r3(xeq,1), r4(xeq,1);];

Cmpc = [1, zeros(1,Lk-1);
    0, 0, 1, zeros(1,Lk-3);
    zeros(1,4), 1, zeros(1,Lk-5);
    zeros(3,8), eye(3), zeros(3,Lk-11);
    zeros(3,14), eye(3), zeros(3,Lk-17);
    lambda_rows;
    zeros(3,Lk-4-3), eye(3), zeros(3,4);
    q2(xeq), zeros(1,Lk-5), q1(xeq), zeros(1,3)];
Cmpc = Sz\Cmpc*Sx;

if xeq(26)<= W.rate_point
    Q_c = [0 1 1 ones(1,3) ones(1,3) 20*ones(1,3) ones(1,3) 5]; % Error Weight (lambda)    
else
    Q_c = [20 1 1 ones(1,3) ones(1,3) zeros(1,3) zeros(1,3) 5]; % Error Weight (omega_r)
end
Qmpc = diag(Q_c);

R_c = 0.1; % Input Weight
Rmpc = R_c*eye(Uk);

R_d = 10; % Change Input Weight
R_del = R_d*eye(Uk);

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
    
    Bi = Ampc*Bi + Bmpc;
    Bcal_u = [Bcal_u; Bi];
end

Bcal_du = Bcal_u;
for i=1:Hu-1
    temp = [zeros(i*Lk,Uk); Bcal_u(1:end-Lk*i,:)];
    Bcal_du = [Bcal_du temp];
end

Ccal = repmat({Cmpc}, 1, Hp);
Ccal = blkdiag(Ccal{:});
Ccal = Ccal(Zk*(Hw-1)+1:end,:); % CHANGE

Psi = Ccal*Acal;
Upsilon = Ccal*Bcal_u;
Theta = Ccal*Bcal_du;

Tau = reshape(refht,[Zk*Hp 1]); % Convert to column vector
Tau = Tau(Zk*(Hw-1)+1:end,:);

Qcal = repmat({Qmpc}, 1, Hp);
Qcal = blkdiag(Qcal{:});
Qcal = Qcal(Zk*(Hw-1)+1:end,:);

Rcal = repmat({Rmpc}, 1, Hu);
Rcal = blkdiag(Rcal{:});

Rcal_del = repmat({R_del}, 1, Hu);
Rcal_del = blkdiag(Rcal_del{:});

ident = eye(Uk);
for i=1:Hu-1
    ident = [ident; eye(Uk)];
end

V = [ident ident];
for i=1:Hu-1
    temp = [zeros(i*Uk,Uk); ident(1:end-Uk*i,:)];
    V = [V temp];
end

%%%%%%%%%%%%
%%% Cost %%%
%%%%%%%%%%%%
deltaU_full = [Uprev; deltaU];
U = V*deltaU_full; % U: u*: The optimal control window
Xcal = Acal*X0 + Bcal_u*Uprev + Bcal_du*deltaU; % P: P*: Optimal states in prediction window
Zcal = Psi*X0 + Upsilon*Uprev + Theta*deltaU;

Epsilon = Tau - Psi*X0 - Upsilon*Uprev;
Gcal = 2*Theta'*Qcal*Epsilon;
Hcal = Theta'*Qcal*Theta + Rcal;

Cost = Epsilon'*Qcal*Epsilon - deltaU'*Gcal + deltaU'*Hcal*deltaU;

% refht_col = reshape(refht,[Zk*Hp 1]); % Convert to column vector
% Cost = Ts*((Zcal-refht_col)'*Qcal*(Zcal-refht_col) + deltaU'*Rcal_del*deltaU +...
%        U'*Rcal*U);

%%%%%%%%%%%%%%%%%%%
%%% Constraints %%%
%%%%%%%%%%%%%%%%%%%

% Constraints = [F*[U;1]<=0; E*[deltaU;1]<=0; G*[Zcal;1]<=0];
Constraints = [F*[U;1]<=0; E*[deltaU;1]<=0];
% Constraints = G*[Zcal;1]<=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLACE YOUR CODE HERE (END)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Yalmip optimizer-object used for simulation and TurtleBot3 control
% ops.mosek.MSK_IPAR_INTPNT_SCALING = 'MSK_SCALING_AGGRESSIVE';
% ops.mosek.MSK_IPAR_SIM_SCALING = 'MSK_SCALING_AGGRESSIVE';
% load('mosekdebug');mosekopt('min write(dump.task.gz)', prob, param)

MPCobj = optimizer(Constraints,Cost,ops,{X0,Uprev,refht},{U,Xcal});
% U: u*: The optimal control window
% P: P*: Optimal states in prediction window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOU CAN PLACE YOUR OWN FUNCTIONS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


