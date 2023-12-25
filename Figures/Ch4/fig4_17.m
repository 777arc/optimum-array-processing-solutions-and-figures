%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.17
% Beam pattern cuts: separable Hamming weighting
% Xin Zhang 9/19/99
% Updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description:
%    Redrawing of Figure 4.17,
% Parameters:
%    N=M=11, dx = dy = lambda/2

close all
clear all

Nx = 11;
Ny = 11;
dx = 0.5;
dy = 0.5;
Dx = [-(Nx-1)/2:1:(Nx-1)/2];
Dy = [-(Ny-1)/2:1:(Ny-1)/2];
ur = -1:1/200:1;

% Hamming
wx = 0.54+0.46*cos(2*pi*Dx.'/Nx);
wy = 0.54+0.46*cos(2*pi*Dy.'/Ny);

phi_set = [0 pi/4 pi/2];

for l = 1:length(phi_set)
   phi = phi_set(l);
   ux = ur*cos(phi);
   uy = ur*sin(phi);
	Ax = exp(j*pi*Dx.'*ux);
	Ay = exp(j*pi*Dy.'*uy);
	Bx = wx'*Ax;
	By = wy'*Ay;
	Beam = Bx.*By;
	BP = abs(Beam)/max(max(abs(Beam)));
   table(l,:) = BP;
end

table = 20*log10(table);
figure;
plot(ur,table(1,:),'b-');
hold on
plot(ur,table(2,:),'g--');
plot(ur,table(3,:),'r-');
hold off
axis([-1 1 -140 0]);
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
h=legend('Cut = 0^o','Cut = 45^o');
set(h,'Fontsize',12)