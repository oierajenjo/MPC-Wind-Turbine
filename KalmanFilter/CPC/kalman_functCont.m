classdef kalman_functCont
    methods
        function res = wp(d)
            res = d(1)*pi/(2*W.L);
        end
        
        function res = cp_ct(la,be,cl,lambdaVec,pitchVec)
            [~,i_la] = min(abs(lambdaVec-abs(la)));
            [~,i_be] = min(abs(pitchVec-be));
            res = cl(i_la,i_be);
        end
        
        function res = ve(x,d)
            res = x(13) + d(1);
        end
        
        function res = Tr(x,d)
            v = ve(x,d);
            res = 0.5*Ae.rho*Ae.Ar*(v-x(3))^3*cp_ct(x(1)*Ae.Rr/(v-x(3)),x(10),cp_l,lambdaVec,pitchVec)/x(1);
        end
        
        function res = Fr(x,d)
            v = ve(x,d);
            res = 0.5*Ae.rho*Ae.Ar*(v-x(3)-x(7))^2*cp_ct(x(1)*Ae.Rr/(v-x(3)),x(10),ct_l,lambdaVec,pitchVec);
        end
        
        function res = f1(x,d)
            tr = Tr(x,d);
            res = (1-D.mu)*tr/(D.Jr+D.Jg) - x(12)/(D.Jr+D.Jg);
        end
        
        function res = f2(x)
            res = x(3);
        end
        
        function res = f3(x)
            res = -(B.B*B.kx + T.k)*x(2)/T.m - (B.B*B.cx + T.c)*x(3)/T.m + B.B*B.kx*x(6)/T.m + B.B*B.cx*x(7)/T.m;
        end
        
        function res = f4(x)
            res = x(5);
        end
        
        function res = f5(x)
            res = 3*x(12)/(2*T.H*T.m) - (B.B*B.ky + T.k)*x(4)/T.m - (B.B*B.cy + T.c)*x(5)/T.m + B.B*B.ky*x(8)/T.m + B.B*B.cy*x(9)/T.m;
        end
        
        function res = f6(x)
            res = x(7);
        end
        
        function res = f7(x,d)
            fr = Fr(x,d);
            res = fr/(B.B*B.m) + B.kx*x(2)/B.m + B.cx*x(3)/B.m - B.kx*x(6)/B.m - B.cx*x(7)/B.m;
        end
        
        function res = f8(x)
            res = x(9);
        end
        
        function res = f9(x)
            res = B.ky*x(4)/B.m + B.cy*x(5)/B.m - B.ky*x(8)/B.m - B.cy*x(9)/B.m;
        end
        
        function res = f10(x)
            res = x(11);
        end
        
        function res = f11(x,u)
            res = Ac.omega^2*u(1) - 2*Ac.omega*Ac.xi*x(11) - Ac.omega^2*x(10);
        end
        
        function res = f12(x,u)
            res = (u(2)-x(12))/Ac.tau;
        end
        
        function res = f13(x,d)
            w = wp(d);
            res = -w*x(13);
        end
        
        function res = f(x,u,d)
            res = x + Ts*[f1(x,d); f2(x); f3(x); f4(x); f5(x); f6(x); f7(x,d);...
                f8(x); f9(x); f10(x); f11(x,u); f12(x,u); f13(x,d)];
        end
        
        function res = h(x,d)
            res = [x(1); f3(x); f5(x); B.l*B.m*f7(x,d); B.l*B.m*f9(x); ...
                D.eta*x(12)*x(1); ve(x,d)-x(3)];
        end
        
        function res = a(d)
            res = 1 - wp(d)*Ts;
        end
        
        function res = sigma_t(d)
            res = ti*d(1)*sqrt((1-a(d)^2)/(1-a(d))^2);
        end
        
        function res = Q(d)
            res = diag([zeros(Lk-1,1); sigma_t(d)^2*wp(d)^2]);
        end
        
        function res = n(d)
            res = sqrt(Q(d))*randn(Lk, 1);
        end
    end
end