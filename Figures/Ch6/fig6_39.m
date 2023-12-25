%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.39
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/30/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Example 6.6.6 , MPDR_DL
%This routine is to test the Average Array Gain, when signal direction 
%is a unform r.v. in [-0.1 0.1]
%The relation between the Loading level and output SNR is studied.

clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';

% uncomment to get figure in text
signalRange = -0.1:1/500:0.1;
%signalRange = -0.1:1/50:0.1;
INRrange = [10 20 30];

SNRrange = [-10:1:30];
LNRrange = -20:50;

ui1 = -0.30;
ui2 = 0.30;

Vi1 = exp(j*n*pi*ui1);
Vi2 = exp(j*n*pi*ui2);

C = exp(j*n*pi*[0]);
f = [1];


k1 = 1;
for INR = 10.^(INRrange/10)
  disp(['loop ' int2str(k1) ' of ' int2str(length(INRrange)) ' ...'])
  k2 = 1;
  for SNR = 10.^(SNRrange/10)
     k3 = 1;
      SINRi = SNR/(1+2*INR);
      for LNR = 10.^(LNRrange/10)
         Gain(k1,k2,k3) = 0;
         for ua = signalRange
      		Vs = exp(j*n*pi*ua);
	   	   Ss = SNR*Vs*Vs';
		      Sn = eye(N) + INR*Vi1*Vi1'+ INR*Vi2*Vi2';
		      Sx = Ss + Sn + LNR*eye(N);     % Diagonal Loading
		      W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %consider White Noise only
		      SINR0 = (W'*Ss*W)/(W'*Sn*W);
		      Gain(k1,k2,k3) = Gain(k1,k2,k3) + SINR0/SINRi;
		   end  % end of ua
		   k3 = k3 + 1;
      end	% end of LNR
   	k2 = k2 + 1;
   end		% end of SNR
	k1 = k1 + 1;
end			% end of INR

Gain = 10*log10(real(Gain)/length(signalRange) ) ;

for k1 = 1:size(Gain,1)
   for k2 = 1:size(Gain,2)
      [A(k1,k2),I(k1,k2)] = max(Gain(k1,k2,:));
   end
end

Load = -21 +  I ;


for h1 = 1:3
   [x,y] = min ( abs(    A(h1,:)-max(A(h1,:)) +3 ));
   c(h1) = -10 + (y-1);
   d(h1) = A(h1, y);
end

plot(SNRrange,A(1,:),'-',SNRrange,A(2,:),'--',SNRrange,A(3,:),'-.');
hold on 
plot(c,d,'-o')
xlabel('{\itSNR} (dB)','Fontsize',14)
ylabel('Optimal gain (dB)','Fontsize',14)
h=legend('{\itINR}=10 dB','{\itINR}=10 dB','{\itINR}=10 dB','3 dB points',3);
set(h,'Fontsize',12)
%title('MPDR\_DL, ui=+/-0.30, ua~U(-0.1, 0.1), N=10')
axis([-10 30 0 50])
grid
    
gtext('{\itLNR}=7 dB','Fontsize',12)
gtext('{\itLNR}=11 dB','Fontsize',12)
gtext('{\itLNR}=15 dB','Fontsize',12)
gtext('{\itLNR}=50 dB','Fontsize',12)
%gtext('LNR=50dB')

gtext('{\itLNR}=11 dB','Fontsize',12)
gtext('{\itLNR}=14 dB','Fontsize',12)
gtext('{\itLNR}=17 dB','Fontsize',12)
gtext('{\itLNR}=20 dB','Fontsize',12)
%gtext('LNR=22dB')

gtext('{\itLNR}=14 dB','Fontsize',12)
gtext('{\itLNR}=17 dB','Fontsize',12)
gtext('{\itLNR}=21 dB','Fontsize',12)
gtext('{\itLNR}=23 dB','Fontsize',12)
%gtext('LNR=26dB')



%save('MPDR')
