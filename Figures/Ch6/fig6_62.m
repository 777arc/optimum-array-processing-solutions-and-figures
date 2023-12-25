%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.62
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
signalRange = -ud:ud/50:ud;
SNRrange = -10:2:40;
LNRrange = -20:2:50;
eigNum = 3;
%SNR = 10^(20/10);
%INR = 10^(30/10);
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
C  = U(:, 1:eigNum);
f = inv(S(1:eigNum,1:eigNum))*C'*P;


k1 = 1;
for INR = 10.^([10 20 30]/10)
   disp(['loop ' int2str(k1) ' of 3 ...'])
   k2 = 1;
   for SNR = 10.^(SNRrange/10)
      k3 = 1;
      for LNR = 10.^(LNRrange/10)
         Gain(k1,k2,k3) = 0;
			for ua = signalRange
			    Vs = exp(j*n*pi*ua);
			    Sn = eye(N) + INR*Vi1*Vi1' + INR*Vi2*Vi2';
			    Ss = SNR*Vs*Vs';
			    Sx = Ss + Sn + LNR*eye(N);
			    W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;
			    SINRo = (W'*Ss*W)/(W'*Sn*W);
			    SINRi = SNR/(2*INR+1);
			    Gain(k1,k2,k3) = Gain(k1,k2,k3) + SINRo/SINRi;
          end		%end of ua
          k3 = k3 + 1;
       end			%end of LNR
       k2 = k2 + 1;
    end				%end of SNR
    k1 = k1 + 1;
end					%end of INR
 
 
Gain = 10*log10( Gain/length(signalRange) );

for k1 = 1:3
   for k2 = 1:length(SNRrange)
      [A(k1,k2),I(k1,k2)] = max( real(Gain(k1,k2,:)) );
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

grid
axis([-10 40 0 50])
xlabel('{\itSNR} (dB)','Fontsize',14)
ylabel('Optimal {\itE}(gain) (dB)','Fontsize',14)
h=legend('{\itINR}=10 dB','{\itINR}=20 dB','{\itINR}=30 dB','3 dB points',3);
set(h,'Fontsize',12)
%title('LCMP\_DL, SNR=20dB, u_i=+/-0.30, 3 Eigen Constraints, B_d(u)=B_c(u), u_d=0.1')
axis([-10 40 0 50])

