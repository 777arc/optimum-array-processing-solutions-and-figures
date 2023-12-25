%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.4
% Last updated 9/5/03 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N=10;
D=[-(N-1)/2:1:(N-1)/2].';
w = [ones(N,1)]/N;

psi = [-1:0.01:1]*pi;
v = exp(j*D*psi);

B = w'*v;

subplot(2,1,1)
alpha = 0;
c = 0.9+0.1*alpha;
plot(psi/pi, 10*log10(abs(B*c).^2))
axis([-1 1 -40 0])
grid on
ylabel('Beampattern(dB)')
title('Problem 2.4.4, N=10, alpha = 0')

subplot(2,1,2)
alpha = 0.9;
c = 0.9+0.1*alpha;
plot(psi/pi, 10*log10(abs(B*c).^2))
axis([-1 1 -40 0])
grid on
ylabel('Beampattern(dB)')
title('alpha = 0.9')

xlabel('y/p','FontName','symbol')

set(gcf,'Paperposition',[0.25 1 8 9])