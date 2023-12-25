%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2.4.5
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

n = [0 1; 0 3; 0 5; 0 8; 1 2; 1 4; 1 8; 2 4; 2 7; 3 4; 3 7; 4 5];

for k = 1:size(n,1)
   B(k,:) = w'*vv-1/N*exp(j*(n(k,1)-M)*psi)-1/N*exp(j*(n(k,2)-M)*psi);
end
B = 10*log10(abs(B).^2);

for l = 1:12
   k = mod(l,3);
   if k == 0
      k = k+3;
   end
   if k == 1
      figure;
   end
   subplot(3,1,k)
   plot(psi/pi,B(l,:))
   axis([-1 1 -40 0])
   grid on
   if k == 3
      xlabel('u')
   end
   ylabel('Beampattern (dB)')
   title(['n1 = ',num2str(n(l,1)),', n2 = ',num2str(n(l,2))])
   if k == 1
   	title(['Problem 2.4.5, N = 10, n1 = ',num2str(n(l,1)),', n2 = ',num2str(n(l,2))])
   end
end
