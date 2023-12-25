%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.40a
%MPDR_DL, with array perturbation
%Relation between signal mitch and array gain is studied
% Derik Lu 
% Lillian Xu updated 04/06/2001
% updated by K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

sigmaP = 0.02;

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
signalRange = 0:BWNN/2/200:BWNN/2;
LNRrange = [10 20 30];
%SNRrange = [-10:5:30];
% uncomment to get figure oin text
SNRrange = [-10:30];
totalTrialNumber = 500;

Ui1 = -0.30;
Ui2 = 0.30;
Us = 0;

C = exp(j*n*pi*[0]);								 %MPDR
f = [1];

INR1 = 10^(30/10);
INR2 = 10^(30/10);

Gain = zeros(length(LNRrange),length(SNRrange));
Gain0 = zeros(length(LNRrange),length(SNRrange));

m = 1;
for LNR = 10.^(LNRrange/10);
   disp(['loop ' int2str(m) ' of 3 ...'])
   k = 1;
   for SNR = 10.^(SNRrange/10)

      for trial = 1:totalTrialNumber
         
         d1 = (n/2+sigmaP*randn(N,1))*Ui1 + sigmaP*randn(N,1)*sqrt(1-Ui1^2);
         d2 = (n/2+sigmaP*randn(N,1))*Ui2 + sigmaP*randn(N,1)*sqrt(1-Ui2^2);
         ds = (n/2+sigmaP*randn(N,1))*Us  + sigmaP*randn(N,1)*sqrt(1-Us^2);
         
         Vi1 = exp(j*2*pi*d1);
         Vi2 = exp(j*2*pi*d2);
         Vs  = exp(j*2*pi*ds);
         Ss = SNR*Vs*Vs';
         Sn = eye(N) + INR1*Vi1*Vi1'+ INR2*Vi2*Vi2';
         Sx = Ss + Sn + LNR*eye(N);     % Diagonal Loading

         W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %consider White Noise only
         
         SINR0 = (W'*Ss*W)/(W'*Sn*W);
         
         SINRi = SNR/(1+INR1+INR2);
         Gain(m,k) = Gain(m,k) + SINR0/SINRi;
      end
      k = k+1;
   end
   m = m + 1;
end


Vs = exp(j*n*pi*Us);
Vi1 = exp(j*n*pi*Ui1);
Vi2 = exp(j*n*pi*Ui2);

m = 1;
for LNR = 10.^(LNRrange/10)
   k = 1;
   for SNR = 10.^(SNRrange/10)
         Ss = SNR*Vs*Vs';
         Sn = eye(N) + INR1*Vi1*Vi1'+ INR2*Vi2*Vi2';
         Sx = Ss + Sn + LNR*eye(N);     % Diagonal Loading
         W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %consider White Noise only
         SINR0 = (W'*Ss*W)/(W'*Sn*W);
         SINRi = SNR/(1+INR1+INR2);
         Gain0(m,k) = SINR0/SINRi;
         k = k+1;
   end
   m = m + 1;
end

Gain = Gain/totalTrialNumber;
Gain = 10*log10(abs(Gain));

Gain0 = 10*log10(abs(Gain0));

x = SNRrange;
plot(x, Gain0(1,:),'-',x,Gain(1,:),'-*')
hold on
plot(x, Gain0(2,:),'--',x,Gain(2,:),'--*')
hold on
plot(x, Gain0(3,:),'-.',x,Gain(3,:),'-.*')
axis([-10 30 0 50])
grid
xlabel('{\it SNR} (dB)','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
%title('MPDR\_DL with perturbation, ua=0, ui=+/-0.30, INR=30dB(each), \sigma_P=0.1*\lambda, 500 trials')
h=legend('{\itLNR}=10 dB, \sigma_{\itp}/\lambda=0','{\itLNR}=10 dB, \sigma_{\itp}/\lambda=0.02', ...
       '{\itLNR}=20 dB, \sigma_{\itp}/\lambda=0','{\itLNR}=20 dB, \sigma_{\itp}/\lambda=0.02', ...
       '{\itLNR}=30 dB, \sigma_{\itp}/\lambda=0','{\itLNR}=30 dB, \sigma_{\itp}/\lambda=0.02',3);
   set(h,'Fontsize',12)






