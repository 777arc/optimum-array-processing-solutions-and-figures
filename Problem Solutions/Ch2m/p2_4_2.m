%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.2
% Xin Zhang
% Last updated 1/8/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 20;
M = (N-1)/2;
n = (-M:M).';
psi = pi*(-1:0.001:1);
vv = exp(j*n*psi);

w = 1/N*[-1*ones(N/2,1); ones(N/2,1)];
B = w'*vv;

subplot(3,1,1)
plot(psi/pi, imag(B))
axis([-1 1 -1.1 1.1])
grid on
ylabel('Im(B)')
title(['Problem 2.4.2, N = ' int2str(N) ', d = \lambda/2'])

subplot(3,1,2)
plot(psi/pi, 10*log10(abs(B).^2))
axis([-1 1 -40 0])
grid on
ylabel('Magnitude (dB)')

subplot(3,1,3)
plot([-1 0 0 1],[-90 -90 90 90])
axis([-1 1 -180 180])
set(gca,'YTick',[-180 -90 0 90 180])
grid on
ylabel('Phase (deg)')
xlabel('\psi/\pi')

