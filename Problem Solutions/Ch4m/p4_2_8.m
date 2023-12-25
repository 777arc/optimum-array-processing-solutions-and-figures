%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 4.2.8
% K. Bell 10/19/99
% updated by K. Bell 10/10/03
% Functions called: bpsphcut_theta, bpsphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

R = 4;
d = 0.4;
N = ceil(pi/asin(0.5*d/R));
d = 2*R*sin(pi/N);
M = floor(2*pi*R);
Np = 2*M+1;

whamm = 0.54+0.46*cos(2*pi*[-M:1:M]'/Np);
theta = pi/2;

h = zeros(Np,1);
for m=-M:1:M
   h(m+M+1) = 1/(((-j)^(m))* besselj(m,2*pi*R*sin(theta)));
end

B1 = exp(-j*2*pi*[0:1:N-1]'*[-M:1:M]/N)/sqrt(N);

w = B1*(h.*whamm);

phi0 = 0;

D = R*[cos(2*pi*[0:1:N-1]'/N) sin(2*pi*[0:1:N-1]'/N) zeros(N,1)];
v0 = exp(j*2*pi*D*[sin(theta)*cos(phi0); sin(theta)*sin(phi0); cos(phi0)]);
w = w/(v0'*w);

bpsphcut_theta(D,w,-60,[90 60 30]*pi/180);
subplot(3,1,1)
title(['Problem 4.2.8, \theta = 90 deg.'])

set(gcf,'paperposition',[0.25 1 8 9])

bpsphere(D,w,-60,100);
title('Problem 4.2.8')
view(60,10)