%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.25
% Magnitude of beam pattern versus ux and uy: rectangular grid array 20*20 with square
%                                             boundary, psi_r=0.4*pi, Kaiser Window
% Xin Zhang 9/19/99	
% updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 20;
beta = 5;
psi_r = 0.4*pi;
[ux,uy] = meshgrid(-1:1/30:1);
psix = pi*ux;
psiy = pi*uy;

Beam = zeros(61,61);
for n = 1:N/2
   for m = 1:N/2
      s1 = n-0.5;
      s2 = m-0.5;
      t1 = sqrt(s1^2+s2^2);
      w = psi_r*besselj(1,psi_r*t1)/2/pi/t1;
      w = w*besseli(0,beta*sqrt(1-(t1/14)^2))/besseli(0,beta);
      B1 = cos(s1*pi*ux);
      B2 = cos(s2*pi*uy);
      Beam = Beam + 4*w*B1.*B2;
   end
end
 
BP = abs(Beam)/max(max(abs(Beam)));

figure
mesh(ux,uy,BP);
grid on
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern','Fontsize',14);
      
      
 
