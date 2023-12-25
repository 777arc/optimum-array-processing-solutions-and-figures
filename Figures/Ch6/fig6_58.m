%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.58
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/16/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 6.7.9  LCMV, No Loading
%This routine is to test the Array Gain, when there is signal mismatch
% Derivative Constraints


clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';
SNR = 10^(20/10);
signalRange = 0:BWNN/2/200:BWNN/2;

C = [ones(N,1),j*n,-n.^2];
m = 1;

for coef = [0 0.4 0.9 1.0]
    f = [1;0;coef*(1-N^2)/12];
    k = 1;
    for ua = signalRange
        Vs = exp(j*n*pi*ua);
        Sn = eye(N);
        Ss = SNR*Vs*Vs';
        W = inv(Sn)*C*inv(C'*inv(Sn)*C)*f;           %consider White Noise only
        SINRo = (W'*Ss*W)/(W'*Sn*W);
        SINRi = SNR;
        A(m,k) = SINRo/SINRi; 								%(6.438)
        k = k+1;
    end
    m = m + 1;
end



A = 10*log10(abs(A));

plot(signalRange/BWNN,A(1,:),'-',signalRange/BWNN,A(2,:),'--',signalRange/BWNN,A(3,:),'-.',signalRange/BWNN,A(4,:),':');
grid
xlabel('{\itu}_{\ita} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
h=legend('{\itg}=0','{\itg}=0.4*{\itB}_{\itc}(0)\prime\prime','{\itg}=0.9*{\itB}_{\itc}(0)\prime\prime','{\itg}={\itB}_{\itc}(0)\prime\prime',3);
set(h,'Fontsize',12)
%title('LCMV, White Noise, 2nd derivative constaints, SNR=20dB')   
axis([0 0.25 4 10])      
