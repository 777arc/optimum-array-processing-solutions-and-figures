%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.45
% Taylor distribution < g(p) > for circular aperture
% Xiaomin Lu 11-2-98	
% Updated by K. Bell 10/2/00	
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
% Data file: bessl.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%

%Parameters:
%   nBar=6, SLL=20dB, scale=4   ;Note: scale=R/lambda


clear all
close all

load('bessl.dat');
Bessel_null = bessl/pi; %data in the table is Bessel_null(u)=0
                         %in the textbook require: Bessel_null(pi*u)=0    
nBar = 4;
SLL = 30;
scale = 5;

A = 1/pi*acosh(10^(SLL/20));
rhu = 0:1/100:1;
p = pi*rhu;
%p=0:pi/100:pi;
weight = 1 ;
for m = 1:nBar-1
      F = -besselj(0,pi*Bessel_null(m));
      for n = 1:nBar-1
         
         Zn = Bessel_null(nBar)*sqrt((A^2+(n-0.5)^2)./(A^2+(nBar-0.5)^2));
         F = F*(1-Bessel_null(m)^2/Zn^2);
         if (n ~= m) F = F/(1-Bessel_null(m)^2/Bessel_null(n)^2); end
      end
      weight = weight+F/besselj(0,pi*Bessel_null(m))^2*besselj(0,Bessel_null(m)*p);

end

 weight = weight*2/pi^2;

weight = abs(weight)/max(abs(weight));

plot(p/pi,weight,'-')
hold on

clear all
N = 2;
rhu = 0:1/100:1;
taper = (1-rhu.^2).^N;
taper = abs(taper)/max(abs(taper));
plot(rhu,taper,'--')


xlabel('\it r/R','Fontsize',14);
ylabel('Weighting','Fontsize',14)
grid
h=legend('Taylor, {\it n_{bar}}=4, SLL=-30 dB','Taper, {\it N}=2');
set(h,'Fontsize',12)



