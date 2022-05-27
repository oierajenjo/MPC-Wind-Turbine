function [xt,yt] = RK4(f,x,u,h,n,v,Ts)
% [Ac,Ae,M,W,B,D,To] = extract_struct(var);
k_1 = f(x,u);
k_2 = f(x+0.5*Ts*k_1,u);
k_3 = f(x+0.5*Ts*k_2,u);
k_4 = f(x+Ts*k_3,u);
xt = x + (1/6)*(k_1+2*k_2+2*k_3+k_4)*Ts + Ts*n;
yt = h(xt) + v;
end