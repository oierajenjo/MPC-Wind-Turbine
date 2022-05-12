% Define object
nlobj = nlmpc(Lk,Zk,Uk);
nlobj.Ts = Ts;
nlobj.PredictionHorizon = Hp;
nlobj.ControlHorizon = Hu;

% Define model
nlobj.Model.StateFcn = f;
nlobj.Model.IsContinuousTime = true;
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


