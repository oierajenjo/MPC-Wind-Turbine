Hp = 5;
Hu = 5;

% Define object
nlobj = nlmpc(Lk,Zk,Uk);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = Hp;
nlobj.ControlHorizon = Hu;

% Define model
nlobj.Model.StateFcn = f;
nlobj.Model.IsContinuousTime = false;
nlmsobj.Model.StateJacFcn = @mystatejac;
nlobj.Model.NumberOfParameters = 1;
nlobj.Model.OutputFcn = hz;

% Specify custom cost function
% nlobj.Optimization.CustomCostFcn = @(X,U,e,data) Ts*sum(sum(U(1:Hp,:)));
% nlobj.Optimization.ReplaceStandardCost = true;

% Specify custom constraint function for the controller
% nlobj.Optimization.CustomEqConFcn = @(X,U,data) X(end,:)';

% Specify costraintes for output variables (OV = z) and manipulated variables (MV = u).
nlobj.Weights.OutputVariables = [20 5 5 5*ones(1,3) 5*ones(1,3) zeros(1,3) zeros(1,6) 5]; % (omega)
nlobj.Weights.ManipulatedVariablesRate = [0.1 0.1 0.1 0.1];

nlobj.OV(1).Min = Ac.omega_min;
nlobj.OV(1).Max = Ac.omega_max;
nlobj.OV(2).Min = To.xd_min;
nlobj.OV(2).Max = To.xd_max;
nlobj.OV(3).Min = To.yd_min;
nlobj.OV(3).Max = To.yd_max;
nlobj.OV(4).Min = B.xd_min;
nlobj.OV(4).Max = B.xd_max;
nlobj.OV(5).Min = B.xd_min;
nlobj.OV(5).Max = B.xd_max;
nlobj.OV(6).Min = B.xd_min;
nlobj.OV(6).Max = B.xd_max;
nlobj.OV(7).Min = B.yd_min;
nlobj.OV(7).Max = B.yd_max;
nlobj.OV(8).Min = B.yd_min;
nlobj.OV(8).Max = B.yd_max;
nlobj.OV(9).Min = B.yd_min;
nlobj.OV(9).Max = B.yd_max;
nlobj.OV(10).Min = W.lambda_min;
nlobj.OV(10).Max = W.lambda_max;
nlobj.OV(11).Min = W.lambda_min;
nlobj.OV(11).Max = W.lambda_max;
nlobj.OV(12).Min = W.lambda_min;
nlobj.OV(12).Max = W.lambda_max;
nlobj.OV(13).Min = Ac.pitch_min;
nlobj.OV(13).Max = Ac.pitch_max;
nlobj.OV(14).Min = Ac.pitch_min;
nlobj.OV(14).Max = Ac.pitch_max;
nlobj.OV(15).Min = Ac.pitch_min;
nlobj.OV(15).Max = Ac.pitch_max;
nlobj.OV(16).Min = Ac.pitch_dot_min;
nlobj.OV(16).Max = Ac.pitch_dot_max;
nlobj.OV(17).Min = Ac.pitch_dot_min;
nlobj.OV(17).Max = Ac.pitch_dot_max;
nlobj.OV(18).Min = Ac.pitch_dot_min;
nlobj.OV(18).Max = Ac.pitch_dot_max;
nlobj.OV(19).Min = Ac.Pe_min;
nlobj.OV(19).Max = Ac.Pe_max;


nlobj.MV(1).Min = Ac.pitch_min;
nlobj.MV(1).Max = Ac.pitch_max;
nlobj.MV(2).Min = Ac.pitch_min;
nlobj.MV(2).Max = Ac.pitch_max;
nlobj.MV(3).Min = Ac.pitch_min;
nlobj.MV(3).Max = Ac.pitch_max;
nlobj.MV(4).Min = Ac.Tg_min;
nlobj.MV(4).Max = Ac.Tg_max;


% Validate prediction model
x0 = x_i;
u0 = u_b(:,1);
validateFcns(nlobj, x0, u0, [], {Ts});

% Create an nlmpcmoveopt object, and specify the sample time parameter.
nloptions = nlmpcmoveopt;
nloptions.Parameters = {Ts};

function [Ampc,Bmpc] = mystatejac(x,~)
A1 = [zeros(1,11), -a1*ones(1,3);
    zeros(1,2), 1, zeros(1,11);
    0, -b3, -b4, zeros(1,2), b1*ones(1,3), b2*ones(1,3), zeros(1,3);
    zeros(1,4), 1, zeros(1,9);
    zeros(1,3), -c4, -c5, zeros(1,6), c2*ones(1,3);
    zeros(3,8), eye(3), zeros(3,3);
    d10(x,0), d3, d4(x,0), zeros(1,2), -d1, zeros(1,2), d2(x,0), zeros(1,5);
    d10(x,1), d3, d4(x,1), zeros(1,3), -d1, zeros(1,2), d2(x,1), zeros(1,4);
    d10(x,2), d3, d4(x,2), zeros(1,4), -d1, zeros(1,2), d2(x,2), zeros(1,3);
    zeros(3,14)];

A2 = [zeros(1,9), -a2, zeros(1,3);
    zeros(3,13);
    c3*ones(1,3), zeros(1,6), -c1, zeros(1,3);
    zeros(3,13);
    d8(x,0), zeros(1,2), d7(x,0), zeros(1,6), d6(x,0), d5(x,0), d9(x,0);
    0, d8(x,1), zeros(1,2), d7(x,1), zeros(1,5), d6(x,1), d5(x,1), d9(x,1);
    0, 0, d8(x,2), zeros(1,2), d7(x,2), zeros(1,4), d6(x,2), d5(x,2), d9(x,2);
    eye(3), zeros(3,10)];

A3 = [-e9(x,0), 0, -e7(x,0), e3, e4, zeros(1,3), -e8(x,0), 0, 0, -e1, 0, 0;
    -e9(x,1), 0, -e7(x,1), e3, e4, zeros(1,4), -e8(x,1), 0, 0, -e1, 0;
    -e9(x,2), 0, -e7(x,2), e3, e4, zeros(1,5), -e8(x,2), 0, 0, -e1;
    zeros(9,14);
    1, zeros(1,13)];

A4 = [-e2(x,0), zeros(1,2), -e11(x,0), zeros(1,6), -e6(x,0), -e5(x,0), -e10(x,0);
    0, -e2(x,1), zeros(1,2), -e11(x,1), zeros(1,5), -e6(x,1), -e5(x,1), -e10(x,1);
    0, 0, -e2(x,2), zeros(1,2), -e11(x,2), zeros(1,4), -e6(x,2), -e5(x,2), -e10(x,2);
    zeros(3,6), eye(3), zeros(3,4);
    zeros(3), -f1*eye(3), -f2*eye(3), zeros(3,4);
    zeros(1,9), -g1, zeros(1,3);
    zeros(1,10), -W.w_p, zeros(1,2);
    zeros(2,13)];

Ampc = eye(Lk)+ Ts*[A1 A2; A3 A4];
%Ampc = Sx\Ampc*Sx;% State Matrix


Bmpc = Ts*[zeros(3,Lk-7), f1*eye(3), zeros(3,4);
    zeros(1,Lk-4), g1, zeros(1,3)]';
%Bmpc = Sx\Bmpc; % Input Matrix
end
