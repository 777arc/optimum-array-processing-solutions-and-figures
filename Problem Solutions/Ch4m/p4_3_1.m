%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 4.3.1
% K. Bell 11/3/99
% updated by K. Bell 11/3/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

u = [-1:0.005:1];
u2 = [0:0.0001:1];
nu = length(u);
nu2 = length(u2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Radial Taper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rlambda = 5;
n=1;
B=ones(1,nu);
B2=ones(1,nu2);
c = ((0.5)^(n+1))/gamma(n+2);

i=find(u);
B(i) = (1/c)*besselj(n+1,2*pi*Rlambda*u(i))./((2*pi*Rlambda*u(i)).^(n+1));

i=find(u2);
B2(i) = (1/c)*besselj(n+1,2*pi*Rlambda*u2(i))./((2*pi*Rlambda*u2(i)).^(n+1));

h1=plot(u, 10*log10(abs(B).^2),'-');

axis([-1 1 -80 0])
ylabel('Beampattern(dB)')
xlabel('u')
title(['Prob. 4.3.1, R=' num2str(Rlambda) '\lambda'])

hold off
grid on

g = find(B2<1/sqrt(2));
HPBW  = 2*u2(min(g));
HPBWn = HPBW*Rlambda*180/pi

g = find(B2<0);
null = u2(min(g));
BWNN = 2*null;
BWNNn = BWNN*Rlambda*180/pi

[b,g] = min(B2);
side = u2(g);
lobe = 10*log10(b^2)

%numerical integral - rectangle rule
% integrate from 0 to pi/2, then double
delta = 0.0001;
theta = [delta/2:delta:0.5]*pi;
uu = sin(theta);
B = (1/c)*besselj(n+1,2*pi*Rlambda*uu)./((2*pi*Rlambda*uu).^(n+1));
G = 0.5*uu.*abs(B).^2;
Dinv = 2*(sum(G)*delta*pi);

D = (1/Dinv)/(2*pi*Rlambda).^2



