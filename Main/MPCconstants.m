%%%%%%%%%%%%%%%
% CONSTRAINTS %
%%%%%%%%%%%%%%%
Sx = eye(Lk);
% Sx(end-3,end-3) = S_means(1); %Tg

Sz = eye(Zk);
% Sz(end,end) = S_means(2); %Pe

% Actuator Range Constraints
fc = [-1; 1];
f_rep = repmat({fc}, 1, Uk*Hu);
F_L = blkdiag(f_rep{:});
u_cons = [Ac.pitch_min; -Ac.pitch_max; Ac.pitch_min; -Ac.pitch_max;
    Ac.pitch_min; -Ac.pitch_max; Ac.Tg_min; -Ac.Tg_max];
U_L = u_cons;
for i=1:Hu-1
    U_L = [U_L; u_cons];
end
F = [F_L U_L];

% Actuator Slew Rate Constraints
ec = [-1; 1];
e_rep = repmat({ec}, 1, Uk-1);
e_rep{end+1} = zeros(2,1);
e_L = blkdiag(e_rep{:});

e_L = repmat({e_L}, 1, Hu);
E_L = blkdiag(e_L{:});
du_cons = [Ac.pitchd_min; -Ac.pitchd_max; Ac.pitchd_min; ...
    -Ac.pitchd_max; Ac.pitchd_min; -Ac.pitchd_max; zeros(2,1)];
dU_L = du_cons;
for i=1:Hu-1
    dU_L = [dU_L; du_cons];
end
E = [E_L dU_L];

% Constraints in Actuator Variables
g = [-1;1];
g_rep = repmat({g}, 1, Zk);
% for i=0:7
%     g_rep{2+i} = zeros(2,1);
% end
% for i=0:2
%     g_rep{14+i} = zeros(2,1);
% end
g_rep{10} = zeros(2,1);
% g_rep{end} = zeros(2,1);
g_L = blkdiag(g_rep{:});

g_L = repmat({g_L}, 1, Hu);
G_L = blkdiag(g_L{:});
z_cons = cell2mat(struct2cell(Z_c));
Z_L = z_cons;
for i=1:Hu-1
    Z_L = [Z_L; z_cons];
end
G = [G_L Z_L];

ops = sdpsettings('solver','mosek','showprogress',0,'verbose',1,...
    'cachesolvers',1,'savedebug',0,'debug',1,'savesolverinput',0,'savesolveroutput',0);

%Define the MPC object with Yalmip:
refht = sdpvar(Zk,Hp);      % Z_{ref}: The window containing the pos-reference
X0 = sdpvar(Lk,1);          % X(k):    The current state
Z0 = sdpvar(Zk,1);          % X(k):    The current state
Uprev = sdpvar(Uk,1);       % U(k-1):  The previous input command.
deltaU = sdpvar(Uk*Hu,1);   % DeltaU

%%%%%%%%%%%%%%%
%%% Lifting %%%
%%%%%%%%%%%%%%%
% Unchanging matrixes
delta_rs = [(Ac.pitch_max-Ac.pitch_min)*ones(1,3) (Ac.Tg_max-Ac.Tg_min)];

% 3 lambdas
% qs = [(Ac.omega_opt-Ac.omega_min) (To.xd_max-To.xd_min) (To.yd_max-To.yd_min)...
%     (B.xd_max-B.xd_min)*ones(1,3) (B.yd_max-B.yd_min)*ones(1,3) ...
%     (W.lambda_max-W.lambda_min)*ones(1,3) (Ac.pitch_max-Ac.pitch_min)*ones(1,3)...
%     (Ac.pitchd_max-Ac.pitchd_min)*ones(1,3) (Ac.Pe_opt-Ac.Pe_min)]; % 3 lambda + xd&yd
% qs = [(Ac.omega_opt-Ac.omega_min) (To.x_max-To.xd_min) (To.y_max-To.yd_min)...
%     (B.x_max-B.x_min)*ones(1,3) (B.y_max-B.y_min)*ones(1,3) ...
%     (W.lambda_max-W.lambda_min)*ones(1,3) (Ac.pitch_max-Ac.pitch_min)*ones(1,3)...
%     (Ac.pitchd_max-Ac.pitchd_min)*ones(1,3) (Ac.Pe_opt-Ac.Pe_min)]; % 3 lambda + x&y

% 1 lambdas
qs = [(Ac.omega_opt-Ac.omega_min) (To.xd_max-To.xd_min) (To.yd_max-To.yd_min)...
    (B.xd_max-B.xd_min)*ones(1,3) (B.yd_max-B.yd_min)*ones(1,3) ...
    (W.lambda_max-W.lambda_min) (Ac.pitch_max-Ac.pitch_min)*ones(1,3)...
    (Ac.pitchd_max-Ac.pitchd_min)*ones(1,3) (Ac.Pe_opt-Ac.Pe_min)]; % 1 lambda + xd&yd
% qs = [(Ac.omega_opt-Ac.omega_min) (To.x_max-To.x_min) (To.y_max-To.y_min)...
%     (B.x_max-B.x_min)*ones(1,3) (B.y_max-B.y_min)*ones(1,3) ...
%     (W.lambda_max-W.lambda_min) (Ac.pitch_max-Ac.pitch_min)*ones(1,3)...
%     (Ac.pitchd_max-Ac.pitchd_min)*ones(1,3) (Ac.Pe_opt-Ac.Pe_min)]; % 1 lambda + x&y

R_c = 0.1; % Input Weight
Rmpc = R_c*eye(Uk);
Rmpc = Rmpc./delta_rs;

R_d = 10; % Change Input Weight
R_del = R_d*eye(Uk);

Rcal = repmat({Rmpc}, 1, Hu);
Rcal = blkdiag(Rcal{:});

Rcal_del = repmat({R_del}, 1, Hu);
Rcal_del = blkdiag(Rcal_del{:});

Tau = reshape(refht,[Zk*Hp 1]); % Convert to column vector
Tau = Tau(Zk*(Hw-1)+1:end,:);

ident = eye(Uk);
for i=1:Hu-1
    ident = [ident; eye(Uk)];
end

V = [ident ident];
for i=1:Hu-1
    temp = [zeros(i*Uk,Uk); ident(1:end-Uk*i,:)];
    V = [V temp];
end

deltaU_full = [Uprev; deltaU];
U = V*deltaU_full; % U: u*: The optimal control window
