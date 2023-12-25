%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.44
% Taylor Beampattern for circular apture
% Xiaomin Lu 11-2-98	
% Updated by K. Bell 10/2/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 9/30/01
% Data file: bessl.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%

%Parameters:
%   nBar=6, SLL=20dB, scale=4   ;Note: scale=R/lambda

clear all
close all

load('bessl.dat');
Bessel_null = bessl/pi; %data in the table is J1(u)=0
                         %in the textbook require: J1(pi*u)=0    
nBar = 6;
SLL = 20;
scale = 5;

u = 0.0000001:1/101:2*scale+2;  % u = 2*R/lambda*sin(theta)

A = 1/pi*acosh(10^(SLL/20));

beam = besselj(1,pi*u)./(pi*u);
for n = 1:nBar-1
   Zn = Bessel_null(nBar)*sqrt((A^2+(n-0.5)^2)./(A^2+(nBar-0.5)^2));
   beam = beam.*(1-u.^2/Zn^2)./(1-u.^2/Bessel_null(n)^2);
end

beam = abs(beam)/max(abs(beam));

plot(u,20*log10(beam))
%title('Taylor pattern of Circ. aperture, nBar=6, SLL=-20dB, R=5*\lambda');
xlabel('\it u_R','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
axis([0 10 -50 0]);
grid


