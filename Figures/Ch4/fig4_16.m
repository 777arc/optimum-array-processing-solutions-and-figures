%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.16
% Beam pattern: separable Hamming weighting
% Xin Zhang 9/18/99
% Updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description:
%    Redrawing of Figure 4.16,
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
ux = [-1:0.025:1];
uy = [-1:0.025:1];
nx = length(ux);
ny = length(uy);

% Hamming
wx = 0.54+0.46*cos(2*pi*Dx.'/Nx);
wy = 0.54+0.46*cos(2*pi*Dy.'/Ny);
Ax = exp(j*pi*Dx.'*ux);
Ay = exp(j*pi*Dy.'*uy);
Bx = wx'*Ax;
By = wy'*Ay;
Bx = ones(nx,1)*Bx;
By = By.'*ones(1,ny);
Beam = Bx.*By;
BP = abs(Beam)/max(max(abs(Beam)));

figure
G = 10*log10(abs(BP).^2);
G = max(G,-150);
mesh(ux,uy,G)
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14);
axis([-1 1 -1 1 -150 0])
