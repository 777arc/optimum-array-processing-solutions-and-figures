%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 6.68 and 6.69
% Xiaomin Lu 6/26/00
% updated by K. Bell 12/1/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/29/01, 11/16/01
% function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis of covariance matrix taper

close all
clear all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
SNR = 10^(0/10);
INR = 10^(30/10);

ui = 0.3;
Vi = exp(j*n*pi*ui);

Va = exp(j*n*pi*0);

BWNN = 4/N;

Ss = SNR*Va*Va';
Sx = Ss + INR*Vi*Vi' + eye(N);

delta_uiRange = -0.1:1/1000:0.1;
u = delta_uiRange;

gammaRange = [0:0.01:0.03];

k1 = 1;
for gamma = gammaRange 
   
   for h1 = 1:N
      for h2 = 1:N
         Sx(h1,h2)  = Sx(h1,h2)*sinc((h1-h2)*gamma/pi);   %here, sinc(x) = sin(x)/x;
      end
   end																%Tapered covariance matrix
   
   W = inv(Sx)*Va/(Va'*inv(Sx)*Va);
   b(k1,:) = W'*exp(j*n*pi*(ui+u));
   
   k2 = 1;
   for uia = ui + delta_uiRange;
   
      Via = exp(j*n*pi*uia);
      Sn = INR*Via*Via' + eye(N);
      SINRo = real(W'*Ss*W)/real(W'*Sn*W);
      SINRi = SNR/(1+INR);
      Gain(k1,k2) = SINRo/SINRi;
      
      k2 = k2 + 1;
   end
   
   k1 = k1 + 1;
end

   
x = delta_uiRange/BWNN;
   
Gain = 10*log10(Gain);
b = 20*log10(abs(b));

plot(x, Gain(1,:),'-',x, Gain(2,:),'--',x, Gain(3,:),'-.',x, Gain(4,:),':')
grid
axis([-0.25 0.25 0 45])
h=legend('\gamma=0','\gamma=0.01','\gamma=0.02','\gamma=0.03',4);
set(h,'Fontsize',12)
xlabel('\Delta{\itu}_{\iti} /{\itBW}_{\itNN}','Fontsize',14)   
ylabel('Gain (dB)','Fontsize',14)
%title('CMT, ua=0(0dB), ui(normial)=0.3(30dB), N=10')



figure
plot(x, b(1,:),'-',x, b(2,:),'--',x, b(3,:),'-.',x, b(4,:),':')
grid
axis([-0.25 0.25 -100 0])
h=legend('\gamma=0','\gamma=0.01','\gamma=0.02','\gamma=0.03',4);
set(h,'Fontsize',12)
xlabel('\Delta{\itu}_{\iti} /{\itBW}_{\itNN}','Fontsize',14)   
ylabel('Beam pattern (dB)','Fontsize',14)
%title('CMT, ua=0(0dB), ui(normial)=0.3(30dB), N=10')
