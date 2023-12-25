%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.59
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/16/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 6.7.10  LCMP, with Loading
%This routine is to test the Array Gain, when there is signal mismatch
% Derivative Constraints


clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
signalRange = 0:BWNN/2/500:BWNN/2;

C = [ones(N,1),j*n,-n.^2];
coef = 0.8;
f = [1;0;coef*(1-N^2)/12];

SNR = 10^(20/10);
INR = 10^(30/10);
Vi1 = exp(j*n*pi*0.30);
Vi2 = exp(-j*n*pi*0.30);

m = 1;
for LNR = 10.^([0 10 20 30]/10)
k = 1;
for ua = signalRange
   Vs = exp(j*n*pi*ua);
   Sn = eye(N) + INR*Vi1*Vi1' + INR*Vi2*Vi2' ;
   Ss = SNR*Vs*Vs';
   Sx = Ss + Sn + LNR*eye(N);							%LCMP with loading
   W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %
   SINRo = (W'*Ss*W)/(W'*Sn*W);
   SINRi = SNR/(2*INR+1);
   A(m,k) = SINRo/SINRi;
   k = k+1;
end
m = m + 1;
end



A = 10*log10(abs(A));

plot(signalRange/BWNN,A(1,:),'-');            %BWnn = 4/N
hold on
plot(signalRange/BWNN,A(2,:),'--');            %BWnn = 4/N
hold on
plot(signalRange/BWNN,A(3,:),'-.');            %BWnn = 4/N
hold on
plot(signalRange/BWNN,A(4,:),':');            %BWnn = 4/N


grid
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
h=legend('{\itLNR}=0 dB','{\itLNR}=10 dB','{\itLNR}=20 dB','{\itLNR}=30 dB',3);
set(h,'Fontsize',12)
%title('LCMP\_DL, SNR=20dB ,Ui=+/-0.30, INR=30dB, g\prime = [1, 0, 0.8*B_c(0)\prime\prime]')   
axis([0 0.25 20 45])