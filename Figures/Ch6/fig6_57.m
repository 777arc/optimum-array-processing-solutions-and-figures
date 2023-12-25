%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.57
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/16/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6.7.8
% LCMV with Directional Constraints.
% Relations between Array Gain, INR, SNR and Optimal LNR are studied

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
ud = 0.1;                       %constraint region

INR = 10^(20/10);
SNRrange = [0 10 20 30];

%LNRrange = -20:50;
LNRrange = -20:2:50;
%signalRange = -ud:ud/200:ud;
signalRange = -ud:ud/50:ud;

ui1 = -0.30;							%interferer configuration
ui2 = 0.30;
Vi1 = exp(j*n*pi*ui1);
Vi2 = exp(j*n*pi*ui2);

uc = 0.0866;
C = exp(j*n*pi*[0 -uc uc]);  	%LCMP with directional constraints
f = [1;sin(-(N/2)*pi*uc)./(N*sin(-.5*pi*uc));sin((N/2)*pi*uc)./(N*sin(.5*pi*uc))];

k1 = 1;
for SNR = 10.^(SNRrange/10)
   disp(['loop ' int2str(k1) ' of 4 ...'])
   k2 = 1;
   for LNR = 10.^(LNRrange/10)
      Gain(k1,k2) = 0;
      for ua = signalRange
         Va = exp(j*n*pi*ua);
         Ss = SNR*Va*Va';
         Sn = eye(N) + INR*Vi1*Vi1' + INR*Vi2*Vi2';
         Sx = Ss + Sn + LNR*eye(N);
         
         W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           	%LCMP
         SINR0 = (W'*Ss*W)/(W'*Sn*W);
         SINRi = SNR/(1+INR+INR);
         Gain(k1,k2) = Gain(k1,k2) + SINR0/SINRi;
      end 				%end of ua
      k2 = k2 + 1;
   end						%end of LNR
   k1 = k1 + 1;
end							%end of SNR

A = 10*log10( abs(Gain)/length(signalRange) );

plot(LNRrange,A(1,:),'-',LNRrange, A(2,:),'--',LNRrange, A(3,:),'-.',LNRrange, A(4,:),':');
xlabel('{\itLNR} (dB)','Fontsize',14)
ylabel('{\itE}(gain) (dB)','Fontsize',14)
%title('LCMP\_DL(DIR), u_i=+/-0.30, u_a~U(-0.1, 0.1), u_c=[0, -0.0866, 0.0866], g=[1, Bc, Bc], INR=20dB')
grid
axis([-20 50 0 50])
h=legend('{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB','{\itSNR}=30 dB');
set(h,'Fontsize',12)
