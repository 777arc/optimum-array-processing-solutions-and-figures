%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.12
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
r = 0.01:0.01:1;         % 1-|rho|^2


INR = 10^(20/10);
Ao = (1+INR)*r;
plot(r,10*log10(Ao),'-');
hold on
INR = 10^(0/10);
Ao = (1+INR)*r;
plot(r,10*log10(Ao),'--');
hold on
INR = 10^(-20/10);
Ao = (1+INR)*r;
plot(r,10*log10(Ao),'-.');
hold on

plot([0.5 0.5], [-20 30],':')
plot([0.95 0.95], [-20 30],':')
xlabel('1-|\rho_{{\its}1}|^{2}','Fontsize',14)
ylabel('{\itA}_{\ito} /{\itN} (dB)','Fontsize',14)


h=legend('{\itINR}=20 dB','{\itINR}=0 dB', '{\itINR}=-20 dB',2);
set(h,'Fontsize',12)
hold off