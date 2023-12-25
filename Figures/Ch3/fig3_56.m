%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.56
% Scanning u-space with beam fans:
% (a) Center beam at u = 0
% (b) Center beam at u = 6/N
% Xin Zhang 4/9/99
% last updated 9/5/00 by K. Bell
% Lillian Xu 04/16/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 32;
n = (-(N-1)/2:(N-1)/2).';
u = -1:0.001:1;
vv = exp(j*pi*n*u);
M = 7;

u0 = [0 -6/N];
for l = 1:length(u0)
   
figure;
hold on;
for m = -(M-1)/2:(M-1)/2
   k = m+(M-1)/2+1;
   w(:,k) = 1/N*exp(-j*pi*n*(u0(l)-2*m/N));
   B(k,:) = w(:,k)'*vv;
   plot(u,real(B(k,:)))
end
hold off;
xlabel('\it u','Fontsize', 14)
ylabel('Beam pattern','Fontsize',14)
%title(['32 elements, 7 beams, center beam at u = ' num2str(u0(l)) ])
if l == 2
 %  title('32 elements, 7 beams, center beam at u = 6/N')
end
grid

end
