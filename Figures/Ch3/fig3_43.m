%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.43
% Beam pattern versus psi:
% Difference beam with uniform weighting, N=10
% Xin Zhang 3/29/99
% Last updated by K. Bell 9/22/00, K. Bell 7/23/01, K. Bell 9/5/01, 9/30/01
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 10;
n = [-(N-1)/2:(N-1)/2].';
d = 0.5;
u = -1:0.001:1;
vv = exp(j*d*2*pi*n*u);
psi = pi*u;

% Symmetric beam pattern
w = 1/N*ones(N,1);
B = w'*vv;

% asymmetric beam pattern
wa1 = 1/N*ones(N/2,1);
wa2 = -1/N*ones(N/2,1);
wa = [wa2; wa1];
Ba = imag(wa'*vv);

% asymmetric beam pattern, 
Ba2 = sinc(N/4*u)./sinc(u/2).*sin(N/4*psi);

figure
h1 = plot(psi,real(B),'--');
hold on;
plot(psi,Ba,'-r');
h2 = plot(psi,Ba2,'-');
hold off;
grid on;
xlabel('{\it \psi}','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)
h=legend([h2 h1],'Difference','Sum');
set(h,'Fontsize',12)