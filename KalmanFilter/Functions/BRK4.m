function [xt,yt] = BRK4(f,xt,h,yt,u_b,d_b,N,n,v,Ts)
% Runge-Kutta 4th order method
for k = 1:N-1
    k_1 = f(xt(:,k), u_b(:,k), d_b(:,k));
    k_2 = f(xt(:,k)+0.5*Ts*k_1, u_b(:,k)+0.5*Ts, d_b(:,k)+0.5*Ts);
    k_3 = f(xt(:,k)+0.5*Ts*k_2, u_b(:,k)+0.5*Ts, d_b(:,k)+0.5*Ts);
    k_4 = f(xt(:,k)+Ts*k_3, u_b(:,k)+Ts, d_b(:,k)+Ts);
    xt(:,k+1) = xt(:,k) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n;  % main equation
    
    yt(:,k) = h(xt(:,k), d_b(:,k)) + v(:,k);
end
yt(:,N) = h(xt(:,N), d_b(:,N)) + v(:,N);
end
