%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.38
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/30/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 6.6.6 , MPDR_DL
%This routine is to test the Array Gain, when there is signal mismatch
%The relation between the Loading level and output SNR is studied.

clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
signalRange = 0.25*BWNN*(0:1/50:1);

ui = [-0.30, 0.30];
Vi = exp(j*n*pi*ui);

C = exp(j*n*pi*[0]);										%MPDR
f = [1];

INR = 10^(30/10);
LNRRange = 10.^([-inf, 20, 30]/10);
SNRRange = 10.^([-10:10:30]/10);

h1 = 1;
for LNR = LNRRange
   h2 = 1;
	for SNR = SNRRange
   	SINRi = SNR/(1+2*INR);
	   h3 = 1;
   	for ua = signalRange
      	Vs = exp(j*n*pi*ua);
	      Ss = SNR*Vs*Vs';
   	   Sn = eye(N) + INR*Vi*Vi';
	      Sx = Ss + Sn + LNR*eye(N);     % Diagonal Loading
   	   W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %consider White Noise only
	      SINR0 = (W'*Ss*W)/(W'*Sn*W);
   	   Gain(h1,h2,h3) = SINR0/SINRi;
      	h3 = h3+1;
      end
      h2 = h2 + 1;
   end
   h1 = h1 + 1;
end

Gain = 10*log10(abs(Gain));

G = squeeze(Gain(1,:,:));
u = signalRange/BWNN;
plot(u,G(1,:),'-',u,G(2,:),'--',u,G(3,:),'-.',u,G(4,:),':',u,G(5,:),'-x');            %BWnn = 4/N
grid
%title('MPDR\_DL with DOA mismatch, um=0, ui=+/-0.30, INR=30dB(each), LNR=0')
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
h=legend('{\itSNR}=-10 dB','{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB','{\itSNR}=30 dB',3);
set(h,'Fontsize',12)
axis([0 0.25 -50 50])      
set(gca, 'xtick',0:0.02:0.25)

figure
G = squeeze(Gain(2,:,:));
u = signalRange/BWNN;
plot(u,G(1,:),'-',u,G(2,:),'--',u,G(3,:),'-.',u,G(4,:),':',u,G(5,:),'-x');            %BWnn = 4/N
grid
%title('MPDR\_DL with DOA mismatch, um=0, ui=+/-0.30, INR=30dB(each), LNR=20dB')
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
h=legend('{\itSNR}=-10 dB','{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB','{\itSNR}=30 dB',3);
set(h,'Fontsize',12)
axis([0 0.25 0 50])      
set(gca, 'xtick',0:0.02:0.25)

figure
G = squeeze(Gain(3,:,:));
u = signalRange/BWNN;
plot(u,G(1,:),'-',u,G(2,:),'--',u,G(3,:),'-.',u,G(4,:),':',u,G(5,:),'-x');            %BWnn = 4/N
grid
%title('MPDR\_DL with DOA mismatch, um=0, ui=+/-0.30, INR=30dB(each), LNR=30dB')
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
h=legend('{\itSNR}=-10 dB','{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB','{\itSNR}=30 dB',3);
set(h,'Fontsize',12)
axis([0 0.25 0 50])      
set(gca, 'xtick',0:0.02:0.25)
