%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.36b
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/30/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ex 6.6.4 
%MPDR, with array perturbation
%Relation between signal mismatch, SNR and Array Gain is studied

clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
% uncomment to get figure in text (takes a long time to run)
SNRrange = [-20:1:50];
%SNRrange = [-20:10:50];

totalTrialNumber = 500;

Ui1 = 0.30;
Us = 0;

INR1 = 10^(30/10);

C = exp(j*n*pi*[0]);								 %MPDR
f = [1];


k1 = 1;
for sigmaP = [0, 0.05, 0.1, 0.2]
   disp(['loop ' int2str(k1) ' of 4 ...'])
   k2 = 1;
   for SNR = 10.^(SNRrange/10)
      Gain(k1,k2) = 0;
      for trial = 1:totalTrialNumber
            ds = (n/2+sigmaP*randn(N,1))*Us  + sigmaP*randn(N,1)*sqrt(1-Us^2);
   	      d1 = (n/2+sigmaP*randn(N,1))*Ui1 + sigmaP*randn(N,1)*sqrt(1-Ui1^2);
	         Vi1 = exp(j*2*pi*d1);
	         Vs  = exp(j*2*pi*ds);
	         Ss = SNR*Vs*Vs';
	         Sn = eye(N) + INR1*Vi1*Vi1';
	         Sx = Ss + Sn;    
	         W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %consider White Noise only
	         SINR0 = (W'*Ss*W)/(W'*Sn*W);
	         SINRi = SNR/(1+INR1);
   	      Gain(k1,k2) = Gain(k1,k2) + SINR0/SINRi;
      end
      k2 = k2 + 1;  
   end
   k1 = k1 + 1;   
end
Gain = 10*log10(real(Gain)/totalTrialNumber);

plot(SNRrange,Gain(1,:),'-',SNRrange,Gain(2,:),'--',SNRrange,Gain(3,:),'-.',SNRrange,Gain(4,:),':')
grid
xlabel('{\itSNR} (dB)','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14);
h=legend('\sigma_{\itp} = 0','\sigma_{\itp} = 0.05\lambda','\sigma_{\itp} = 0.10\lambda','\sigma_{\itp} = 0.20\lambda',3);
set(h,'Fontsize',12)
axis([-20 50 -80 60])
%title(['MPDR with perturbation., ui=',num2str(Ui1),...
%      '(',num2str(10*log10(INR1)),'dB), ua=',num2str(Us),...
%      ', N=',num2str(N)])
