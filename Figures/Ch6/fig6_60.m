%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.60
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
signalRange = -0.1:1/100:0.1;
SNRrange = -10:2:40;
INRrange = [10 20 30];
LNRrange = -20:2:50;

C = [ones(N,1),j*n,-n.^2];
coef = 0.8;
f = [1;0;coef*(1-N^2)/12];

Vi1 = exp(j*n*pi*0.30);
Vi2 = exp(-j*n*pi*0.30);

k1 = 1;
for INR = 10.^(INRrange/10)
    disp(['loop ' int2str(k1) ' of 3 ...'])
    k2 = 1;
    for SNR = 10.^(SNRrange/10)
        k3 = 1;
        SINRi = SNR/(2*INR+1);
        for LNR = 10.^(LNRrange/10)
            Gain(k1,k2,k3) = 0;
            for ua = signalRange
                Vs = exp(j*n*pi*ua);
                Sn = eye(N) + INR*Vi1*Vi1' + INR*Vi2*Vi2' ;
                Ss = SNR*Vs*Vs';
                Sx = Ss + Sn + LNR*eye(N);							%LCMP with loading
                W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %
                SINRo = (W'*Ss*W)/(W'*Sn*W);
                Gain(k1,k2,k3) = Gain(k1,k2,k3) + SINRo/SINRi;
            end	%end of ua
            k3 = k3 + 1;
        end		%end of LNR
        k2 = k2 + 1;
    end			%end of SNR
    k1 = k1 + 1;
end				%end of INR


Gain = 10*log10( abs(Gain)/length(signalRange) );   

for k1 = 1:length(INRrange)
   for k2 = 1:length(SNRrange)
      [A(k1,k2), I(k1,k2)] = max( Gain(k1,k2,:) );
   end
end

for h1 = 1:3
   [x,y] = min ( abs(    A(h1,:)-max(A(h1,:)) +3 ));
   c(h1) = -10 + 2*(y-1);
   d(h1) = A(h1, y);
end

plot(SNRrange,A(1,:),'-',SNRrange,A(2,:),'--',SNRrange,A(3,:),'-.')
hold on
plot(c,d,'-o')

axis([-10 40 0 50])
xlabel('{\itSNR} (dB)','Fontsize',14)
ylabel('Optimal {\itE}(gain) (dB)','Fontsize',14)
h=legend('{\itINR}=10 dB','{\itINR}=20 dB','{\itINR}=30 dB','3 dB points',3);
set(h,'Fontsize',12)
%title('LCMP\_DL(2nd DER), ui=+/-0.30, ua~U(-0.1, 0.1), N=10, g = [1; 0; 0.8*B_C(0)\prime\prime]')
grid

