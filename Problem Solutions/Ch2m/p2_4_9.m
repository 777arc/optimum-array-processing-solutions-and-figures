%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.9
% Xin Zhang
% Last updated 1/15/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 8;
M = (N-1)/2;
n = (-M:M).';
d = 5/8;
u = -2:0.001:2;
vv = exp(j*2*pi*d*n*u);
w = 1/N*[ones(N,1)];
B8 = w'*vv;
B8 = 10*log10(abs(B8).^2);

N = 10;
M = (N-1)/2;
n = (-M:M).';
d = 1/2;
vv = exp(j*2*pi*d*n*u);
w = 1/N*[ones(N,1)];
B10 = w'*vv;
B10 = 10*log10(abs(B10).^2);

plot(u,B8,'-');
hold on
plot(u,B10,'--');
axis([-2 2 -25 0])
grid on
xlabel('u')
ylabel('Beampattern (dB)')
title('Problem 2.4.9')
legend('N=8, d=5/8\lambda','N=10, d=1/2\lambda')
hold off
