%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.13
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
r = 0.01:0.01:1;         % 1-|rho|^2


N = 10;
n = (-(N-1)/2:(N-1)/2)';
BWNN = 4/N;
vs = ones(N,1);

ui = BWNN*(0:0.01:2.5);
ni = length(ui);

SNR = 1;

for INR = 10.^([0 20]/10)
    vi = exp(j*n*pi*ui);
    rho_sq = abs(vs'*vi/N).^2;
    Ao = N*((1+INR)/(1+N*INR))*(ones(1,ni)+N*INR*(ones(1,ni)-rho_sq));
    Ac = N*(1+INR)*(1+N*INR*rho_sq).^(-1);
    h1=plot(ui/BWNN,10*log10(Ao),'-');
    hold on
    h2=plot(ui/BWNN,10*log10(Ao./Ac),'--');
end

plot([0.23 0.23], [0 35],':')
plot([0.5 0.5], [0 35],':')
xlabel('{\itu}_{\itI} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('dB','Fontsize',14)

h=legend([h1 h2],'{\itA}_{\ito}','{\itA}_{\ito} /{\itA}_{\itc}',0);
set(h,'Fontsize',12)
axis([0 2.5 0 35])
hold off