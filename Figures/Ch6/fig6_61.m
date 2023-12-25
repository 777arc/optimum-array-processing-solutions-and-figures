%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.61
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/16/01
% function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 6.7.11
% LCMP ,eigenvector constraints, with diagonal loading.

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
ud = 0.1;                       %constraint region
BWNN = 4/N;
signalRange = 0:BWNN/2/200:BWNN/2;
eigNum = 3;

SNR = 10^(20/10);
INR = 10^(30/10);
Vi1 = exp(j*n*0.3*pi);
Vi2 = exp(-j*n*0.3*pi);

t = -ud:1/500:ud;
H = exp(j*n*pi*t)*diag(sinc(-(N/2)*t)./(sinc(-.5*t)));
H = transpose(H);
P = transpose(1/500*trapz(H));				%Ud = Bc

for m = 1:N
   for k = 1:N
      Q(m,k) = 2*ud*sinc(ud*(m-k));    
   end
end
[U,S,V] = svd(Q);
C = U(:,1:eigNum);
f = inv(S(1:eigNum,1:eigNum))*C'*P;

m = 1;
for LNR = 10.^([0 10 20 30]/10)
   
   k = 1;
   for ua = signalRange
      Vs = exp(j*n*pi*ua);
      Sn = eye(N) + INR*Vi1*Vi1' + INR*Vi2*Vi2';
      Ss = SNR*Vs*Vs';
      Sx = Ss + Sn + LNR*eye(N);
      W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;
      SINRo = (W'*Ss*W)/(W'*Sn*W);
      SINRi = SNR/(2*INR+1);
      A(m,k) = SINRo/SINRi;
      k = k + 1;
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
h=legend('{\itLNR}=0 dB','{\itLNR}=10 dB','{\itLNR}=20 dB','{\itLNR}=30 dB');
set(h,'Fontsize',12)
%title('LCMP\_DL, SNR=20dB, INR=30dB, U_i=+/-0.30, 3 Eigen Constraints, Bd(u)=Bc(u), ud=0.1')   
axis([0 0.25 20 50])






   

   
   


   
