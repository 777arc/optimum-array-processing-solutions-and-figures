%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.26
%  K. Bell 7/25/01, 9/5/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 11;
n = (-(N-1)/2:(N-1)/2)';
us = -0.7;
vs = exp(j*n*pi*us);

u = -1:0.001:1;
v = exp(j*n*pi*u);
%SNR = 1;
%INR = 10.^([0 10 20]/10);

a1= -0.9;
phi=[0.5 -0.4 -0.6 -0.67 -0.7];
psi = [-1:0.01:1]*pi;
nf = size(psi,2);

n1 = size(phi,2);

for n=1:n1
    z1 = -a1*exp(j*pi*phi(n));
    P = ones(1,nf)./((abs(ones(1,nf)-z1*exp(-j*psi)) ).^2);
    figure
    subplot(2,1,1)
    plot(psi/pi,10*log10(P))
    ylabel('dB','Fontsize',14)
    r = [1+abs(a1)^2 a1*exp(-j*pi*phi(n)) zeros(1,N-2)];
    c = [1+abs(a1)^2 a1*exp(j*pi*phi(n)) zeros(1,N-2)];
    Sxinv = toeplitz(c,r);
    Sxinv(1,1)=1;
    Sxinv(N,N) = 1;
    w = inv(vs'*Sxinv*vs)*Sxinv*vs;
    B = w'*v;
    
    subplot(2,1,2)
    plot(u,10*log10(abs(B).^2),'-')
    hold on
    plot(-0.7*[1 1],[-50 10],'-')
    plot(phi(n)*[1 1],[-50 10],'--')
    axis([-1 1 -45 5])
    xlabel('\psi/\pi','Fontsize',14)
    ylabel('Beam pattern (dB)','Fontsize',14)
    grid on
end
