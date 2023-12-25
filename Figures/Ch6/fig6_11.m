%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.11
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
r = 0:0.01:1;         % 1-|rho|^2

for INR = 10.^([-20 0 20]/10)
    Ao = N*(1+INR)*(1+N*INR*r)/(1+N*INR);
    Ac = N*(1+INR)*(1+N*INR*(1-r)).^(-1);
    h1=plot(r,10*log10(Ac),'-');
    hold on
    h2=plot(r,10*log10(Ao),'--');
end

plot([0.5 0.5], [0 35],':')
plot([0.95 0.95], [0 35],':')
xlabel('1-|\rho_{{\its}1}|^{2}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)

text(0.1,25,'{\itINR}=20 dB','Fontsize',12)
text(0.6,3,'{\itINR}=20 dB','Fontsize',12)
text(0.2,5,'{\itINR}=0 dB','Fontsize',12)
text(0.6,13,'{\itINR}=0 dB','Fontsize',12)
text(0.1,10.5,'{\itINR}=-20 dB','Fontsize',12)

h=legend([h1 h2],'Conventional','Optimum',2);
set(h,'Fontsize',12)
hold off