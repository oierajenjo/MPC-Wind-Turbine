clear all
close all

load('BladedFiles\performancemap_data.mat')

theta = pitchVec;
lambda = lambdaVec;

theta_int = interp1(1:length(theta),theta,1:0.562:length(theta));

[dCp_dTheta, dCp_dLambda] = gradient(cp_l);

theta_eq = 0.2;
lambda_eq = 10;

[theta_rest, theta_index] = min(abs(theta-theta_eq));
[lambda_rest, lambda_index] = min(abs(lambda-lambda_eq));

dCp_dTheta(lambda_index,theta_index)
dCp_dLambda(lambda_index,theta_index)