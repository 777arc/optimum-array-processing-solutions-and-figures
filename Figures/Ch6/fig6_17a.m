%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.17a
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/19/00, 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 10;
alphaRange = 0:1/500:1;
INRRange = 10.^([-10 0 10 20]/10);

h1 = 1;
for INR = INRRange
   h2 = 1;
   for alpha2 = alphaRange
      T(h1,h2) = 1/N*( 1 + 2*N*INR*alpha2 + N^2*INR^2*alpha2)/( 1 + N*INR*alpha2)^2;
      h2 = h2 + 1;
   end
   h1 = h1 + 1;
end

u = alphaRange;
T = 10*log10(T);
plot(u, T(1,:),'-',u, T(2,:),'--',u, T(3,:),'-.', u,T(4,:),':')
xlabel('\alpha^2','FontSize',14)
ylabel('{\itT} (dB)','FontSize',14)
axis([0 1 -20 20])
grid
h=legend('{\itINR}=-10 dB','{\itINR}=0 dB','{\itINR}=10 dB','{\itINR}=20 dB');
set(h,'FontSize',12)




