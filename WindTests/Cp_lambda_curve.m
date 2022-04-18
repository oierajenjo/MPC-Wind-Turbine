clear all
close all

load('BladedFiles\performancemap_data.mat')

i = 1;
for pitch = 0:(5*pi/180):(20*pi/180)
    for lamb = 0:0.01:20
        cp(i) = cp_ct(lamb,pitch,cp_l,lambdaVec,pitchVec);
        i = i+1;
    end
    i = 1;
    figure(1)
    plot(0:0.01:20,cp,'LineWidth',1)
    hold on
end
xlabel('lambda [-]')
ylabel('C_p [-]')
legend('theta=0º','theta=5º','theta=10º','theta=15º','theta=20º')
ylim([-0.5 0.5])
hold off