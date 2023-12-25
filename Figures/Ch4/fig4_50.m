%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.50
% Difference Pattern for a Circular Aperture with 
%  the Generic Aperture Distribution
% Xiaomin Lu 11-2-98	
% Updated by K. Bell 10/2/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%

%   
clear all
close all


u = 0:1/301:7;
h = 0.0001;
t = 1/2*(besselj(0,pi*u)-bessel(2,pi*u)); %differentiation of J1(v) 
D = u.*t./(1-(u/0.5860670).^2);
D = abs(D)/max(abs(D));
plot(u,20*log10(D))
grid
axis([0 7 -50 0])
%title('Beampattern of Circular Aperture with the Generic Distribution')
xlabel('\it u_R','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
