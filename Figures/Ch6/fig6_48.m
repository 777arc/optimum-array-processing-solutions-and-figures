%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.48
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/16/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 6.7.6 LCMP beamformer, No Loading
%This routine is to test the Array Gain, when there is signal mismatch
%

clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
signalRange = 0.25*BWNN*(0:1/100:1);


detU = 0.0866;
C = exp(j*n*pi*[0 -detU detU]);
f = [1;sin(-(N/2)*pi*detU)./(N*sin(-.5*pi*detU));sin((N/2)*pi*detU)./(N*sin(.5*pi*detU))];

m = 1;
for SNR = 10.^([-20:10:30]/10)
   k = 1;
   for ua = signalRange
      Vs = exp(j*n*pi*ua);
      Ss = SNR*Vs*Vs';
      Sn = eye(N);
      Sx = Ss + Sn;
      W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %LCMP, white noise only
      A(m,k) = ((W'*Ss*W)/(W'*Sn*W) )/ SNR; % white noise
      k = k+1;
   end
   m = m + 1;
end

A = 10*log10(abs(A));
plot(signalRange/BWNN,A(1,:),'-')
hold on
plot(signalRange/BWNN,A(2,:),'--')
hold on
plot(signalRange/BWNN,A(3,:),'-.')
hold on
plot(signalRange/BWNN,A(4,:),':')
hold on
plot(signalRange/BWNN,A(5,:),'-x')
hold on
plot(signalRange/BWNN,A(6,:),'-')

grid
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
h=legend('{\itSNR}=-20 dB','{\itSNR}=-10 dB','{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB','{\itSNR}=30 dB',3);
set(h,'Fontsize',12)
%title('LCMP, White Noise, uc=[0 +/-0.0866], g=[1;Bc;Bc]')
axis([0 0.25 -50 20])      
line([0.0866/BWNN 0.0866/BWNN],[-80 15])   
text(0.0866/BWNN+0.01,14,'{\itu}_{{\itc}1}','Fontsize',12)
set(gca,'xtick',0:0.02:0.25)
   
