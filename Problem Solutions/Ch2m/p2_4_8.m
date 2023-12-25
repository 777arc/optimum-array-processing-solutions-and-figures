%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.8
% Xin Zhang
% Last updated 1/25/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = (2:100);
HPBW = 0.886*2./N;
u4 = sqrt(3./(N.^2-1));

figure;
plot(N,HPBW,'-',N,u4,'--')
title('Problem 2.4.8, d = \lambda/2')
xlabel('N')
legend('HPBW','\Deltau4')
