%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.12
% Polar plot of uniformly weighted rectangular array
% Xiaomin Lu  11/2/98	
% Updated by K. Bell 9/29/00
% Function called: sinc, polardb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters:
%    N=M=10, phi = 0, pi/4, pi/2


close all
clear all
N = 10;
M = 10;

theta = -pi:pi/200:pi;
i = 1;
for phi = [0, pi/4];
   %for phi = pi/4
   ux = sin(theta).*cos(phi);
   uy = sin(theta).*sin(phi);
   Beam = sinc(N*ux/2)./sinc(ux/2);
   Beam = Beam.*sinc(M*uy/2)./sinc(uy/2);
   Beam = abs(Beam)/max(max(abs(Beam)));
   Beam = 20*log10(Beam);
   table(i,:) = Beam;
   i = i+1;
end


polardb(theta,table(1,:),-50,'-');
%text(-6,6.4,'Pattern cut of uniform planar array, N=M=10, dx=dy=lambda/2, phi=0');
text(-0.3,-6.5,'(a)');
figure
polardb(theta,table(2,:),-80,'-');
%text(-10,10,'Pattern cut of uniform planar array, N=M=10, dx=dy=lambda/2, phi=pi/4');
text(-0.3,-10.3,'(b)');
