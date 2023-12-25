%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.6
% Last updated 9/5/03 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N=7;
D=[-(N-1)/2:1:(N-1)/2].';
w = [ones(N,1)]/N;

psi = [-1:0.01:1]*pi;
v = exp(j*D*psi);

B = w'*v;

D4 = [-3 -2 1 3]';
w4 = ones(4,1)/4;
v4 = exp(j*D4*psi);

B4 = w4'*v4;

h2=plot(psi/pi, 10*log10(abs(B).^2),'--');
hold on
h1 =plot(psi/pi, 10*log10(abs(B4).^2),'-');

axis([-1 1 -40 0])
grid on
ylabel('Beampattern(dB)')
title('Problem 2.4.6')
legend([h1 h2],'4 elements','7 elements')
xlabel('y/p','FontName','symbol')

hold off

u = [0.2:0.0001:0.3];
vv = exp(j*D4*pi*u);
B = w4'*vv;
[y,I]=min(abs(B));
umin = u(I)
BWmm_u = 2*umin
psimin = pi*umin
BWmm_psi = 2*psimin
Bmin = y

