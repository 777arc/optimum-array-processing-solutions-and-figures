%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.63
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/29/01, 11/16/01
% function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following methods are considered:
%
%	 (1)  Orthogonal Quiescent Pattern(Q: conventional beampattern)
%
%
%	 System configuration:
%			10 element linear array

clear all
close all

N = 10;							%Array size	

n = (-(N-1)/2:(N-1)/2)';
ud = 0.1;
w0 = ones(N,1)/N;				%conventional beampattern pointed at broadside(u=0)
Ne = 2;							%Number of principle eigenvectors chosed

signalRange = -ud:0.02:ud;
SNRrange = -10:40;
INRrange = 10:10:30;
LNRrange = -20:50;


ui1 = -0.30;
ui2 = 0.30;

Vi1 = exp(j*n*pi*ui1);
Vi2 = exp(j*n*pi*ui2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% QP2 constraints with conventional beampattern

wdBar = w0/(w0'*w0);		
PcWdBar = eye(N) - wdBar*inv(wdBar'*wdBar)*wdBar';
for h1 = 1:N
	for h2 = 1:N
   	Rs(h1,h2) = 2*ud*sinc( (h1-h2)*ud );
   end
end
RsBar = PcWdBar*Rs*PcWdBar;
[U,lambda,V] = svd(RsBar);
Cs = V(:,1:Ne);					%select 3 eigenvectors
C = [wdBar,Cs];
f = [1;zeros(Ne,1)];

%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1 = 1;
for INR = 10.^(INRrange/10)
   disp(['loop ' int2str(k1) ' of 3 ...'])
   k2 = 1;
   for SNR = 10.^(SNRrange/10)
      k3 = 1;
      for LNR = 10.^(LNRrange/10)
         Gain(k1,k2,k3) = 0;
         for ua = signalRange
            Va = exp(j*n*pi*ua);
            Ss = SNR*Va*Va';
            Sn = eye(N) + INR*Vi1*Vi1' + INR*Vi2*Vi2';
            Sx = Ss + Sn + LNR*eye(N);
 			   W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           %
			   SINRo = (W'*Ss*W)/(W'*Sn*W);
			   SINRi = SNR/(2*INR+1);
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

Load = -20 + 1*(I-1);

for h1 = 1:3
   [x,y] = min ( abs(    A(h1,:)-max(A(h1,:)) +3 ));
   c(h1) = -10 + 1*(y-1);
   d(h1) = A(h1, y);
end

%save fig6_63
plot(SNRrange,A(1,:),'-',SNRrange,A(2,:),'--',SNRrange,A(3,:),'-.')
hold on
plot(c,d,'-o')
axis([-10 40 0 50])
xlabel('{\itSNR} (dB)','Fontsize',14)
ylabel('Optimal {\itE}(gain) (dB)','Fontsize',14)
h=legend('{\itINR}=10 dB','{\itINR}=20 dB','{\itINR}=30 dB','3 dB points',3);
set(h,'Fontsize',12)
%title('QP2, ui=+/-0.30, ua~U(-0.1, 0.1), N=10, Ne=2, Quiescent: Conv')
grid
x = [11 21 31 41];
plot(-11+x,A(1,x),'x')
plot(-11+x,A(2,x),'x')
plot(-11+x,A(3,x),'x')
gtext('{\itLNR}=-6 dB','Fontsize',12)
gtext('{\itLNR}=3 dB','Fontsize',12)
gtext('{\itLNR}=9 dB','Fontsize',12)
gtext('{\itLNR}=32 dB','Fontsize',12)

gtext('{\itLNR}=-1 dB','Fontsize',12)
gtext('{\itLNR}=6 dB','Fontsize',12)
gtext('{\itLNR}=11 dB','Fontsize',12)
gtext('{\itLNR}=15 dB','Fontsize',12)

gtext('{\itLNR}=3 dB','Fontsize',12)
gtext('{\itLNR}=9 dB','Fontsize',12)
gtext('{\itLNR}=14 dB','Fontsize',12)
gtext('{\itLNR}=18 dB','Fontsize',12)
hold off
      

        
        
