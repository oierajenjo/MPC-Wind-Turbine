function [xt,yt] = RK4(f,xt,u_b,h,yt,n,v,var,Ts)
disp('Running True Values')
% [Ac,Ae,M,W,B,D,To] = extract_struct(var);

N = size(u_b,2);
% Runge-Kutta 4th order method
for k = 1:N-1
 % main equation
    k_1 = f(xt(:,k),u_b(:,k));
    k_2 = f(xt(:,k)+0.5*Ts*k_1,u_b(:,k)+0.5*Ts);
    k_3 = f(xt(:,k)+0.5*Ts*k_2,u_b(:,k)+0.5*Ts);
    k_4 = f(xt(:,k)+Ts*k_3,u_b(:,k)+Ts);
    xt(:,k+1) = xt(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n(xt(:,k)); 
    yt(:,k) = h(xt(:,k)) + v(:,k);
end
yt(:,N) = h(xt(:,N)) + v(:,N);
end