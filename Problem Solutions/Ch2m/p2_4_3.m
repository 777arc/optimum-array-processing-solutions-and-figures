%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.3
% Xin Zhang
% Last updated 2/11/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 10;
M = (N-1)/2;
D = (-M:M).';
psi = pi*(-1:0.001:1);
vv = exp(j*D*psi);

w = 1/N*[ones(N,1)];

for n = 0:5
   B(n+1,:) = w'*vv-1/N*exp(j*(n-M)*psi);
end
B = 10*log10(abs(B).^2);

for l = 1:6
   if l > 3
      k = l-3;
   else
      k = l;
   end
   if l == 1 | l == 4
      figure;
   end
   subplot(3,1,k)
   plot(psi/pi,B(l,:))
   axis([-1 1 -40 0])
   grid on
   if l == 3 | l == 6
      xlabel('u')
   end
   ylabel('Beampattern (dB)')
   title(['n = ',num2str(l-1)])
   if l == 1 | l == 4
   	title(['Problem 2.4.3, N = 10, n = ',num2str(l-1)'])
   end
end
