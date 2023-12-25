%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.46
% Comparative Beampatterns of Taylor and Taper weighting, R=5*\lambda
% Xiaomin Lu  11-2-98	
% Updated by K. Bell 10/2/00 
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
% Data file: bessl.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%


clear all
close all

load('bessl.dat');
Bessel_null = bessl/pi; %data in the table is J1(u)=0
                         %in the textbook require: J1(pi*u)=0    

nBar = 4;
SLL = 30;
scale = 5;
u = 0.00001:1/100:2*scale;  % u = 2*R/lambda*sin(theta)

A = 1/pi*acosh(10^(SLL/20));


beam = besselj(1,pi*u)./(pi*u);
for n = 1:nBar-1
   Zn = Bessel_null(nBar)*sqrt((A^2+(n-0.5)^2)./(A^2+(nBar-0.5)^2));
   beam = beam.*(1-u.^2/Zn^2)./(1-u.^2/Bessel_null(n)^2);
end

beam = abs(beam)/max(abs(beam));
%subplot(2,1,1);

plot(u,20*log10(beam))
hold on

clear all
N = 2;
scale = 5;
i = 1;
for u = 0:1/50:2*scale
   p = 0:pi/100:pi;
   y = p.*(1-(p/pi).^2).^N.*besselj(0,u*p);
   beam(i) = trapz(p,y);
   i = i+1;
end
x = 0:1/50:2*scale;
 
beam = abs(beam)/max(abs(beam));

plot(x,20*log10(beam),'--')
%title('Comparative Beampatterns of Taylor and Taper weighting, R=5*\lambda');
axis([0 10 -80 0]);
grid
xlabel('\it u_R','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14)
h=legend('Taylor, {\it n_{bar}}=4, SLL=-30 dB','Taper, {\it N}=2');
set(h,'Fontsize',12)

