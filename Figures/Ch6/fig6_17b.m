%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.17b
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/19/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
BWNN = 4/N;
uiRange = BWNN*(0:1/500:1);
INRRange = 10.^([-10 0 10 20]/10);

SNR = 1;
Va = ones(N,1);

h1 = 1;
for INR = INRRange
   h2 = 1;
   for ui = uiRange
      Vi = exp(j*n*pi*ui);
      Sx = SNR*Va*Va' + INR*Vi*Vi' + eye(N);
      W = inv(Sx)*Va/(Va'*inv(Sx)*Va);
      T(h1,h2) = norm(W)^2;
      h2 = h2 + 1;
   end
   h1 = h1 + 1;
end

u = uiRange/BWNN;
T = 10*log10(T);
plot(u, T(1,:),'-',u, T(2,:),'--',u, T(3,:),'-.', u,T(4,:),':')
xlabel('{\itu}_{\itI} /{\itBW}_{\itNN}','FontSize',14)
ylabel('{\itT} (dB)','FontSize',14)
axis([0 1 -20 20])
grid
h=legend('{\itINR}=-10 dB','{\itINR}=0 dB','{\itINR}=10 dB','{\itINR}=20 dB');
set(h,'FontSize',12)





