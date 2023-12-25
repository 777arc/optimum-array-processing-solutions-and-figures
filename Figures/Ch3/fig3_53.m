%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.53
% Beam patterns: N = 11
% (a) u-space
% (b) theta-space
% Xin Zhang 4/9/99
% last updated 9/5/00 by K. Bell
% Lillian Xu 04/16/2001, 7/23/01, 9/30/01
% function called: polardb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 11;
M = (N-1)/2;
n = (-M:M).';
u = -1:0.001:1;
vv = exp(j*pi*n*u);
u0 = 0;
figure;
hold on;
for m = -M:M
   k = m+M+1;
   w(:,k) = 1/N*exp(-j*pi*n*(u0-2*m/N));
   B(k,:) = w(:,k)'*vv;
   plot(u,real(B(k,:)))
%   pause
end
hold off;
xlabel('\it u','Fontsize', 14)
ylabel('Beam pattern','Fontsize',14)
%title('N = 11')
grid

theta = pi*(-1:0.001:1);
vt = exp(j*pi*n*cos(theta));
lim = -40;
figure;
for m = -M:M
   k = m+M+1;
   Bt(k,:) = w(:,k)'*vt;
   rho(k,:) = 20*log10(abs(Bt(k,:)));
   rp(k,:) = rho(k,:);
   for l = 1:length(theta)
      if (cos(theta(l)) < u0+2*m/N-2/N) | (cos(theta(l)) > u0+2*m/N+2/N)
         rp(k,l) = -100;
      end
   end
   hpol(k,:) = polardb(theta,rp(k,:),lim);
   hold on;
end
hold off;
