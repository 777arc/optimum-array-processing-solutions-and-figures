%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.6.12
% Xin Zhang 2/12/99
% Last updated by K. Bell 9/13/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 10;
M = (N-1)/2;
n = (-M:M).';
u = -1:0.001:1;

w = 1/N*[ones(N,1)];
vv = exp(j*pi*n*u);
B_n = w'*vv;

sigma_p = 0.01;
sls = (2*pi*sigma_p)^2;
p = exp(-sls*(1-u.^2));
BP = (abs(B_n).^2).*p+(1-p)/N;
h1 = plot(u, 10*log10(BP),'-');
hold on

sigma_p = 0.05;
sls = (2*pi*sigma_p)^2;
p = exp(-sls*(1-u.^2));
BP = (abs(B_n).^2).*p+(1-p)/N;
h2 = plot(u, 10*log10(BP),'--');

sigma_p = 0.1;
sls = (2*pi*sigma_p)^2;
p = exp(-sls*(1-u.^2));
BP = (abs(B_n).^2).*p+(1-p)/N;
h3 = plot(u, 10*log10(BP),'-.');

sigma_p = 0.5;
sls = (2*pi*sigma_p)^2;
p = exp(-sls*(1-u.^2));
BP = (abs(B_n).^2).*p+(1-p)/N;
h4 = plot(u, 10*log10(BP),':');

axis([-1 1 -40 0])
ylabel('Beampattern (dB)')
title('Problem 2.6.12, perturbations in y-axis, \phi = \pi/2')
xlabel('u')
legend('\sigma_\lambda = 2\pi(0.01)', '\sigma_\lambda = 2\pi(0.05)', ...
       '\sigma_\lambda = 2\pi(0.1)', '\sigma_\lambda = 2\pi(0.5)');
hold off

%%%%%%%%%%%
figure
sigma_p = 0.1;
sls = (2*pi*sigma_p)^2;

phi = 0/180*pi;
p = exp(-sls*(1-u.^2)*sin(phi)^2);
BP = (abs(B_n).^2).*p+(1-p)/N;
h1 = plot(u, 10*log10(BP),'-');
hold on

phi = 30/180*pi;
p = exp(-sls*(1-u.^2)*sin(phi)^2);
BP = (abs(B_n).^2).*p+(1-p)/N;
h2 = plot(u, 10*log10(BP),'--');

phi = 60/180*pi;
p = exp(-sls*(1-u.^2)*sin(phi)^2);
BP = (abs(B_n).^2).*p+(1-p)/N;
h3 = plot(u, 10*log10(BP),'-.');

phi = 90/180*pi;
p = exp(-sls*(1-u.^2)*sin(phi)^2);
BP = (abs(B_n).^2).*p+(1-p)/N;
h4 = plot(u, 10*log10(BP),':');

axis([-1 1 -40 0])
ylabel('Beampattern (dB)')
title('Problem 2.6.12, perturbations in y-axis, \sigma_\lambda = 2\pi(0.1)')
xlabel('u')
legend('\phi = 0','\phi = \pi/6','\phi = \pi/3','\phi = \pi/2');
hold off
