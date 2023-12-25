%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.8
% Beam patterns for Kaiser weighting
% beta = 3, N = 11, 21, and 41
% Xin Zhang
% Last updated 04/16/2001 Lillian Xu, K. Bell 7/20/01,9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N1 = 11;
N2 = 21;
N3 = 41;
beta = 3;

bes = besseli(0,beta);
for n = -(N1-1)/2:(N1-1)/2
   x = 1-(2*n/N1)^2;
   wn1(n+1+(N1-1)/2) = besseli(0,beta*sqrt(x));
end
wn1 = wn1'/bes;

for n = -(N2-1)/2:(N2-1)/2
   x = 1-(2*n/N2)^2;
   wn2(n+1+(N2-1)/2) = besseli(0,beta*sqrt(x));
end
wn2 = wn2'/bes;

for n = -(N3-1)/2:(N3-1)/2
   x = 1-(2*n/N3)^2;
   wn3(n+1+(N3-1)/2) = besseli(0,beta*sqrt(x));
end
wn3 = wn3'/bes;

u=(-1:0.001:1);
amf1=[-(N1-1)/2:(N1-1)/2]';
Vo1=exp(i*amf1*pi*u);
amf2=[-(N2-1)/2:(N2-1)/2]';
Vo2=exp(i*amf2*pi*u);
amf3=[-(N3-1)/2:(N3-1)/2]';
Vo3=exp(i*amf3*pi*u);
wn1=wn1/sum(wn1);
wn2=wn2/sum(wn2);
wn3=wn3/sum(wn3);
Beam1=20*log10(abs(wn1'*Vo1));
Beam2=20*log10(abs(wn2'*Vo2));
Beam3=20*log10(abs(wn3'*Vo3));
plot(u,Beam3,'-',u,Beam2,'--',u,Beam1,'-.')
h=legend('{\it N} = 41','{\it N} = 21','{\it N} = 11');
axis([-1 1 -50 0])
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
