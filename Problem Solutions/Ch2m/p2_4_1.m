%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.1
% Xin Zhang
% Updated 9/8/99 by K. Bell
% Last updated 9/13/00 by K. Bell
% Functions called: polardb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 21;
M = (N-1)/2;
n = (-M:M).';
ang = pi*(-1:0.001:1);
u = cos(ang);
vv = exp(j*pi*n*u);
theta = [0 15 30 60 90]/180*pi;

for k = 1:length(theta)
   w = 1/N*ones(N,1).*exp(j*n*pi*cos(theta(k)));
   B(k,:) = w'*vv;
end

B = 10*log10(abs(B).^2);

for k = 1:length(theta)
   figure
   h=polardb(ang,B(k,:),-40);
   title(['Problem 2.4.1, N = 21, \theta_T = ' num2str(theta(k)/pi*180) ' (deg)'])
   hold off
end

scan_limit = (acos(1-0.450*2/N))*180/pi

bw = zeros(1,5);
for k = 1:length(theta)
   ur(k) = cos(theta(k))+0.450*2/N;
   if k == 1
      bw(k) = 2*acos(1-0.450*2/N);
   else
      if abs(ur(k)) <= 1
      theta_L = acos(ur(k));
      theta_R = acos(cos(theta(k))-0.450*2/N);
      bw(k) = theta_R - theta_L;
      end
   end
end
% HPBW
bw = bw/pi*180