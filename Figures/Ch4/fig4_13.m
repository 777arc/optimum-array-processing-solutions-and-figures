%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.13
% Cut plot of uniformly weighted rectangular array
% Xiaomin Lu 11/2/98	
% Updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parameters:
%    N=M=10, dx = dy = lambda/2 , betax = betay = 0

close all
clear all
N = 10;
M = 10;

i = 1;
ur = -1:1/100:1;
for phi = [0 pi/4]
   ux = ur*cos(phi);
   uy = ur*sin(phi);
   
   %[ux,uy] = meshgrid(-1:1/50:1);
   psix = pi*ux;
   psiy = pi*uy;

   Beam = sinc(1/pi*N*psix/2)./sinc(1/pi*psix/2);
   Beam = Beam.*sinc(1/pi*M*psiy/2)./sinc(1/pi*psiy/2);
   Beam = abs(Beam)/max(abs(Beam));
   table(i,:) = Beam;
   i = i+1;
end

plot(ur,20*log10(table(1,:)));
axis([-1 1 -50 0]);
grid
%title('Plot cut of Uniform Planar Array, N=M=10, d=\lambda/2, \phi=0')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
text(-0.027,-56,'(a)','Fontsize',14)

figure
plot(ur,20*log10(table(2,:)));
axis([-1 1 -50 0]);
grid
%title('Plot cut of Uniform Planar Array, N=M=10, d=\lambda/2, \phi=\pi/4 ')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
text(-0.027,-56,'(b)','Fontsize',14)



