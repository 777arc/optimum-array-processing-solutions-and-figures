%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.8.4
% Xin Zhang 2/13/99
% Updated 9/13/00 by K. Bell
% Functions called: bpsphere, sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

Nx = 10;
Nz = 10;
nx = (-(Nx-1)/2:(Nx-1)/2).';
nz = (-(Nz-1)/2:(Nz-1)/2).';
dx = 1/2;
dz = 1/2;
theta = pi*(0:0.001:1);
%phi = [0 15 30 45 60 75 90]*pi/180;
phi = [90 75 60 45 30 0]*pi/180;

for l = 1:length(phi)
   
   % Calculate the beam pattern
   for k = 1:length(theta)
    	B_x(l,k) = sinc(Nx/2*sin(theta(k))*cos(phi(l)))/sinc(1/2*sin(theta(k))*cos(phi(l)));
      B_z(l,k) = sinc(Nz/2*cos(theta(k)))/sinc(1/2*cos(theta(k)));
      B(l,k) = B_x(l,k)*B_z(l,k);
      B(l,k) = 20*log10(abs(B(l,k)));
      if B(l,k) < -60
         B(l,k) = -60;
      end
   end
   
   % Output the results
   r = mod(l+1,2)+1;
   if r == 1
      figure;
   end
   subplot(2,1,r)
   plot(cos(theta),B(l,:))
   ylabel('Beampattern (dB)')
   grid on;
   if r == 1
      title(['Problem 2.8.4, \phi = ',num2str(phi(l)/pi*180),' degrees'])
   else
      title(['\phi = ',num2str(phi(l)/pi*180),' degrees'])
	   xlabel('cos(\theta)')
   end
end

% 3D plot
% Uniform Planar Array in x-z plane
Dx = [-(Nx-1)/2:1:(Nx-1)/2].'*dx;
Dz = [-(Nz-1)/2:1:(Nz-1)/2].'*dz;

D=[kron(Dx,ones(Nz,1)) [zeros(Nx*Nz,1)] kron(ones(Nx,1),Dz)];
w = ones(Nx*Nz,1)/(Nx*Nz);
lim = -60;
m=100;

bpsphere(D,w,lim,m)
title('Problem 2.8.4, 10x10 planar array in x-z plane')