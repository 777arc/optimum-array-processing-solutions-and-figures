%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.35
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/30/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
BWNN = 4/N;
signalRange = (0:0.25/100:0.25)*BWNN;

Vm = exp(j*n*pi*0);

k1 = 1;
for ASNR = 0:10:30
   M = 10^(ASNR/10);
   k2 = 1;
   for ua = signalRange
      Va = exp(j*n*pi*ua);
      Bc = 1/N*Vm'*Va;
      A(k1,k2) = N*abs(Bc)^2/( 1 + ( 2*M + M^2 )*( 1 - abs(Bc)^2 )  );
      k2 = k2 + 1;
   end
   k1  = k1 + 1;
end

A = 10*log10(A);
plot(signalRange/BWNN, A(1,:),'-',signalRange/BWNN, A(2,:),'--',signalRange/BWNN, A(3,:),'-.')
hold on
plot(signalRange/BWNN, A(4,:),':')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gain coresponding to CRB
hold on
 
h2 = 1;
for K = [2*N, 10*N]
   h1 = 1;
   for ASNR = 0:2:40
	   M = 10^(ASNR/10);
	   a = 3/(8*K);
   	b = inv( 1 - 1/N^2);
	   delta(h1) = sqrt(a*(1/M + 1/M^2)*b);
   
	%   [x,I] = min( abs(signalRange/BWNN-delta(h1)) )
   	Va = exp( j*n*pi*delta(h1)*BWNN );
	   Bc = 1/N*Vm'*Va;
   	CR_Gain(h2,h1) = N*abs(Bc)^2/( 1 + ( 2*M + M^2 )*( 1 - abs(Bc)^2 )  );
%   CR_Gain(h1) = A(h1,I);
   
	   h1 = h1 + 1;
   end
   h2 = h2 + 1;
   
end

CR_Gain = 10*log10(CR_Gain);
plot(delta(4:length(delta)),CR_Gain(1,4:length(delta)), '-o',delta(6:length(delta)),CR_Gain(2,6:length(delta)), '-p');
axis([0 0.25 -60 10])
grid
h=legend('{\itASNR}=0 dB','{\itASNR}=10 dB','{\itASNR}=20 dB','{\itASNR}=30 dB','CRB,{\itK}=2{\itN}','CRB,{\itK}=10{\itN}',4);
set(h,'Fontsize',12)
xlabel('{\itu}_{\ite} /{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
%title('MPDR with DOA mismatch, N=10')
axis([0 0.1 -40 20])
set(gca, 'xtick',[0:0.01:0.25])
set(gca, 'ytick',[-40:5:20])

   
   




