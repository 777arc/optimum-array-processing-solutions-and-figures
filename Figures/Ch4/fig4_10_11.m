%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.10 & 4.11
% Beampattern of uniformly weighted rectangular array						
% Xiaomin Lu 11/2/98	
% Updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 9/30/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description:
%    Redrawing of Figure 4.10, 
% Parameters:
%    N=M=10, dx = dy = lambda/2 , betax = betay = 0

close all
clear all

N = 10;
M = 10;

[ux,uy] = meshgrid(-1:1/80:1);
psix = pi*ux;
psiy = pi*uy;

Beam = sinc(1/pi*N*psix/2)./sinc(1/pi*psix/2);
Beam = Beam.*sinc(1/pi*M*psiy/2)./sinc(1/pi*psiy/2);
Beam = abs(Beam)/max(max(abs(Beam)));

mesh(ux,uy,Beam);
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern','Fontsize',14);
%title('Response of planar array, N=M=10, dx=dy=\lambda/2, uniform weight')

figure
Beam = 20*log10(Beam);
for i = 1:size(Beam,1)
   for j = 1:size(Beam,2)
      if (Beam(i,j)< -60)
         Beam(i,j) = -60;
      end
   end
end

mesh(ux,uy,Beam)
  
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14);
%title('Reponse of planar array, N=M=10, dx=dy=\lambda/2, uniform weight')



