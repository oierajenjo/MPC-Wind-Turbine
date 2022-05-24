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

Ampc = [A1 A2; A3 A4];
% eig(Ampc)
Ampc = eye(Lk) + Ts*Ampc;
% Ampc = Sx\Ampc*Sx;% State Matrix


Bmpc = Ts*[zeros(3,Lk-7), f1*eye(3), zeros(3,4);
    zeros(1,Lk-4), g1, zeros(1,3)]';
% Bmpc = Sx\Bmpc; % Input Matrix


% lambda_rows = [
%     r1(xeq,0), 0, r6(xeq,0), zeros(1,5), r7(xeq,0), zeros(1,5), r2(xeq,0), zeros(1,9), r5(xeq,0), r3(xeq,0), r4(xeq,0);
%     r1(xeq,1), 0, r6(xeq,1), zeros(1,6), r7(xeq,1), zeros(1,5), r2(xeq,1), zeros(1,8), r5(xeq,1), r3(xeq,1), r4(xeq,1);
%     r1(xeq,2), 0, r6(xeq,2), zeros(1,7), r7(xeq,2), zeros(1,5), r2(xeq,2), zeros(1,7), r5(xeq,2), r3(xeq,2), r4(xeq,2)];

% Cmpc = [1, zeros(1,Lk-1);
%     0, 0, 1, zeros(1,Lk-3);
%     zeros(1,4), 1, zeros(1,Lk-5);
%     zeros(3,8), eye(3), zeros(3,Lk-11);
%     zeros(3,14), eye(3), zeros(3,Lk-17);
%     lambda_rows;
%     zeros(6,Lk-4-6), eye(6), zeros(6,4);
%     q2(xeq), zeros(1,Lk-5), q1(xeq), zeros(1,3)]; % 3 lambdas + xd&yd
% Cmpc = Sz\Cmpc*Sx;

% Cmpc = [1, zeros(1,Lk-1);
%     0, 1, 0, zeros(1,Lk-3);
%     zeros(1,3), 1, zeros(1,Lk-4);
%     zeros(3,5), eye(3), zeros(3,Lk-8);
%     zeros(3,11), eye(3), zeros(3,Lk-14);
%     lambda_rows;
%     zeros(6,Lk-4-6), eye(6), zeros(6,4);
%     q2(xeq), zeros(1,Lk-5), q1(xeq), zeros(1,3)]; % 3 lambdas + x&y
% Cmpc = Sz\Cmpc*Sx;


lambda_row = [rr1(xeq), 0, rr4(xeq), zeros(1,Lk-6), rr3(xeq), rr2(xeq), 0];

Cmpc = [1, zeros(1,Lk-1);
    0, 0, 1, zeros(1,Lk-3);
    zeros(1,4), 1, zeros(1,Lk-5);
    zeros(3,8), eye(3), zeros(3,Lk-11);
    zeros(3,14), eye(3), zeros(3,Lk-17);
    lambda_row;
    zeros(6,Lk-4-6), eye(6), zeros(6,4);
    q2(xeq), zeros(1,Lk-5), q1(xeq), zeros(1,3)]; % 1 lambda + xd&yd

% Cmpc = [1, zeros(1,Lk-1);
%     0, 1, 0, zeros(1,Lk-3);
%     zeros(1,3), 1, zeros(1,Lk-4);
%     zeros(3,5), eye(3), zeros(3,Lk-8);
%     zeros(3,11), eye(3), zeros(3,Lk-14);
%     lambda_row;
%     zeros(6,Lk-4-6), eye(6), zeros(6,4);
%     q2(xeq), zeros(1,Lk-5), q1(xeq), zeros(1,3)]; % 1 lambda + x&y


% Cy = [1, zeros(1,Lk-1);
%     0, -b3, -b4, zeros(1,2), b1*ones(1,3), b2*ones(1,3), zeros(1,Lk-11);
%     zeros(1,3), -c4, -c5, zeros(1,6), c2*ones(1,3), c3*ones(1,3), zeros(1,6), -c1, zeros(1,3);
%     zeros(3,5), p1*eye(3), zeros(3,Lk-8);
%     zeros(3,11), p2*eye(3), zeros(3,Lk-14);
%     zeros(3,17), eye(3), zeros(3,Lk-20);
%     q2(xeq), zeros(1,Lk-5), q1(xeq), zeros(1,3);
%     zeros(1,2), -1, zeros(1,Lk-6), 1, 1, 0;
%     zeros(1,Lk-1), 1];


% if xeq(26)<= W.rate_point % 3 lambda
%     Q_c = [0 10 10 10*ones(1,3) 10*ones(1,3) 20*ones(1,3) 20*ones(1,3) zeros(1,3) 0]./qs; % Error Weight (lambda)
% else
%     Q_c = [20 10 10 10*ones(1,3) 10*ones(1,3) zeros(1,3) zeros(1,3) zeros(1,3) 20]./qs; % Error Weight (omega_r)
% end
if xeq(Lk-1)<= W.rate_point % 1 lambdas
    Q_c = [10 5 5 5*ones(1,3) 5*ones(1,3) 20 20*ones(1,3) zeros(1,3) 0]./qs; % Error Weight (lambda)
else
    Q_c = [20 5 5 5*ones(1,3) 5*ones(1,3) 0 zeros(1,3) zeros(1,3) 20]./qs; % Error Weight (omega_r)
end
Qmpc = diag(Q_c);


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
Ccal = Ccal(Zk*(Hw-1)+1:end,:);

Qcal = repmat({Qmpc}, 1, Hp);
Qcal = blkdiag(Qcal{:});
Qcal = Qcal(Zk*(Hw-1)+1:end,:);


%%%%%%%%%%%%
%%% Cost %%%
%%%%%%%%%%%%
Psi = Ccal*Acal;
Upsilon = Ccal*Bcal_u;
Theta = Ccal*Bcal_du;

Xcal = Acal*X0 + Bcal_u*Uprev + Bcal_du*deltaU; % P: P*: Optimal states in prediction window
Zcal = Psi*X0 + Upsilon*Uprev + Theta*deltaU;

Epsilon = Tau - Psi*X0 - Upsilon*Uprev;
Gcal = 2*Theta'*Qcal*Epsilon;
Hcal = Theta'*Qcal*Theta + Rcal;

Cost = Ts*(Epsilon'*Qcal*Epsilon - deltaU'*Gcal + deltaU'*Hcal*deltaU);

% refht_col = reshape(refht,[Zk*Hp 1]); % Convert to column vector
% Cost = Ts*((Zcal-refht_col)'*Qcal*(Zcal-refht_col) + deltaU'*Rcal_del*deltaU +...
%        U'*Rcal*U);


%%%%%%%%%%%%%%%%%%%
%%% Constraints %%%
%%%%%%%%%%%%%%%%%%%

% Constraints = [F*[U;1]<=0; G*[Zcal;1]<=0; E*[deltaU;1]<=0];
Constraints = [F*[U;1]<=0; E*[deltaU;1]<=0];
% Constraints = [F*[U;1]<=0; G*[Zcal;1]<=0];
% Constraints = G*[Zcal;1]<=0;

% The Yalmip optimizer-object used for simulation and control
MPCobj = optimizer(Constraints,Cost,ops,{X0,Uprev,refht},{U,Xcal,Zcal});
% U: u*: The optimal control window
% P: P*: Optimal states in prediction window
