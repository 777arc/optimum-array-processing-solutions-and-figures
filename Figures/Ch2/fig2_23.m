%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.23
% Beam pattern for 10-element uniform array (d=lambda/2)
% scanned to 30 degrees (60 degrees from broadside)
% Xin Zhang 1/20/99
% Last updated by K. Bell 6/25/01
% Functions called: polardb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2).';
theta = pi*(-1:0.001:1);
u = cos(theta);
d = 1/2;
vv = exp(j*2*pi*d*n*u);
theta_T = 30/180*pi;
w = 1/N*ones(N,1).*exp(j*n*pi*cos(theta_T));
B = w'*vv;
B = 10*log10(abs(B).^2);
figure
h=polardb(theta,B,-40);
hold off
