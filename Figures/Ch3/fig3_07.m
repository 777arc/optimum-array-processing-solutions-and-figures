%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.7
% Kaiser weighting
% (a)weighting for beta = 3 and 6
% (b)Beam patterns for beta = 3 and 6
% Xin Zhang updated 3/17/99
% Lillian Xiaolan Xu updated 09/20/2000
% Lillian Xiaolan Xu updated 04/16/2001, K. Bell 7/20/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 11;
Zl = (N-1)/2;
beta1 = 3;
beta2 = 6;

for n = -(N-1)/2:(N-1)/2
   x = 1-(2*n/N)^2;
   wn1(n+1+(N-1)/2,1) = besseli(0,beta1*sqrt(x));
end

for n = -(N-1)/2:(N-1)/2
   x = 1-(2*n/N)^2;
   wn2(n+1+(N-1)/2,1) = besseli(0,beta2*sqrt(x));
end

nn = -Zl/2:1/2:Zl/2;

u=(-1:0.001:1);
amf=[-(N-1)/2:(N-1)/2]';
Vo=exp(i*amf*pi*u);
w1=wn1/sum(wn1);
w2=wn2/sum(wn2);
Beam1=20*log10(abs(w1'*Vo));
Beam2=20*log10(abs(w2'*Vo));

bes1 = besseli(0,beta1);
bes2 = besseli(0,beta2);
wn1 = wn1/bes1;
wn2 = wn2/bes2;

figure
subplot(2,1,1)
plot(nn,wn1,'-',nn,wn2,'--')
xlabel('{\it z}/\lambda','Fontsize',14)
ylabel('Weights','Fontsize',14)
%title('Kaiser Window for 11 elements')
h=legend('\beta = 3','\beta = 6');
set(h,'Fontsize',12)
hold on;
plot(nn,wn1,'o')
plot(nn,wn2,'o')
axis([-2.5 2.5 -0.02 1.02])
hold off;
text(-0.0625,-0.28,'(a)','Fontsize',14)
subplot(2,1,2)
plot(u,Beam1,'-',u,Beam2,'--')
h=legend('\beta = 3','\beta = 6');
axis([-1 1 -60 0])
text(-0.025,-75,'(b)','Fontsize',14)
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
