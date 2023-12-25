%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.54
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/16/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 6.7.8 LCMV beamformer
%This routine is to test the Array Gain, when there is signal mismatch


clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
SNRrange = 0:10:30;
signalRange = 0:BWNN/2/200:BWNN/2;

detU = 0.0866;
C = exp(j*n*pi*[0 -detU detU]);
f = [1;sin(-(N/2)*pi*detU)./(N*sin(-.5*pi*detU));sin((N/2)*pi*detU)./(N*sin(.5*pi*detU))];
Vi1 = exp(j*n*pi*0.30);
Vi2 = exp(-j*n*pi*0.30);
INR = 10^(10/10);
LNR = 10^(15/10);

k1 = 1;
for SNR = 10.^(SNRrange/10)
    k2 = 1;
    for ua = 0:BWNN/2/200:BWNN/2
        Va = exp(j*n*pi*ua);				
        Ss = SNR*Va*Va';
        Sn = INR*Vi1*Vi1' + INR*Vi2*Vi2' + eye(N);	      %LCMP
        Sx = Ss + Sn + LNR*eye(N);
        W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %
        SINRo = (W'*Ss*W)/(W'*Sn*W); 
        SINRi = SNR/(1+2*INR);
        A(k1,k2) = SINRo/SINRi;
        k2 = k2+1;
    end		
    k1 = k1 + 1;
end


A = 10*log10(abs(A));
x = signalRange/BWNN;
plot(x,A(1,:),'-',x,A(2,:),'--',x,A(3,:),'-.',x,A(4,:),':');
%BWnn = 4/N
grid
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
h=legend('{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB','{\itSNR}=30 dB',3);
set(h,'Fontsize',12)
%title('LCMP\_DL, ui=+/-0.30(10dB each), LNR=10dB, uc=[0,+/-0.0866], g=[1;Bc;Bc]')
set(gca,'xtick',0:0.02:0.25)
axis([0 0.25 -10 30])      
line([0.0866/BWNN 0.0866/BWNN],[-50 45])   
text(0.0866/BWNN+0.002,25,'{\itu}_{{\itc}1}','Fontsize',12)
   
