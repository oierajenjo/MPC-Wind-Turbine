function [xt,yt] = RK4_NL(f,x,u,h,n,v,Ts)
% [Ac,Ae,M,W,B,D,To] = extract_struct(var);
k_1 = f(x,u,Ts);
k_2 = f(x+0.5*Ts*k_1,u+0.5*Ts,Ts);
k_3 = f(x+0.5*Ts*k_2,u+0.5*Ts,Ts);
k_4 = f(x+Ts*k_3,u+Ts,Ts);
xt = x + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n;
yt = h(xt) + v;
end