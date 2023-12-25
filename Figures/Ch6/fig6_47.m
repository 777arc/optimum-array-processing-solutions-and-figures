%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.47
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/16/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 6.7.5 LCMV beamformer
%This routine is to test the Array Gain, when there is signal mismatch


clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
SNR = 10^(20/10);
SINRi = SNR;

signalRange = 0.25*BWNN*(0:1/100:1);

detU = 0.0866;
C = exp(j*n*pi*[0 -detU detU]);
f1 = [1;1;1];
f2 = [1;sin(-(N/2)*pi*detU)./(N*sin(-.5*pi*detU));sin((N/2)*pi*detU)./(N*sin(.5*pi*detU))];
f3= [1;1/2+1/2*sin(-(N/2)*pi*detU)./(N*sin(-.5*pi*detU));1/2+1/2*sin((N/2)*pi*detU)./(N*sin(.5*pi*detU))];
f4= [1;2/3+1/3*sin(-(N/2)*pi*detU)./(N*sin(-.5*pi*detU));2/3+1/3*sin((N/2)*pi*detU)./(N*sin(.5*pi*detU))];
f = [f1,f2,f3,f4];

for k1 = 1:4
    g = f(:,k1);
    k2 = 1;
    for ua = signalRange
        Vs = exp(j*n*pi*ua);				
        Sn = eye(N);											%LCMV, here only white noise is available
        W = inv(Sn)*C*inv(C'*inv(Sn)*C)*g;           %
        Ss = SNR*Vs*Vs';
        A(k1,k2) = ((W'*Ss*W)/(W'*Sn*W))/SINRi; % white noise
        k2 = k2+1;
    end		
end


A = 10*log10(abs(A));
x = signalRange/BWNN;
plot(x,A(1,:),'-',x,A(2,:),'--',x,A(3,:),'-.',x,A(4,:),':');            %BWnn = 4/N
grid
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
%title('LCMV-DIR,WhiteNoise, SNR=20dB, uc=[0,+/-0.0866]')
axis([0 0.25 4 10])
h=legend('{\itg}_{{\itc}1}','{\itg}_{{\itc}2}','{\itg}_{{\itc}3}','{\itg}_{{\itc}4}',3);
set(h,'Fontsize',12)
line([0.0866/BWNN 0.0866/BWNN],[0 9])   
text(0.0866/BWNN+0.01,9,'{\itu}_{{\itc}1}','Fontsize',12)
set(gca,'xtick',0:0.02:0.25)
set(gca, 'ytick',4:0.5:10)
