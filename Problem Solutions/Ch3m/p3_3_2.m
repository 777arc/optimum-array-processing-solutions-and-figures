%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3.3.2
% K. Bell 10/13/98
% Updated 9/22/00 by K. Bell
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=21;
d=0.5;
D = d*[-(N-1)/2:1:(N-1)/2].';

u = [-1:0.001:1];
v = exp(j*2*pi*D*u);

Bd = 0.342*((sqrt(1-u.^2)).^(-1)).*((u>=0.5)&u<=cos(pi*20/180));

subplot(2,1,1)
plot(u,Bd)
hold on
delta = 0.3;

uk = [-(N-1)/2:1:(N-1)/2]*2/N+delta*2/N;
Bp = 0.342*((sqrt(1-uk.^2)).^(-1)).*((uk>=0.5)&uk<=cos(pi*20/180));
plot(uk,Bp,'o')

Bc = zeros(1,length(u));

for k=1:N
   Bc = Bc+Bp(k)*sinc((u-uk(k))*N/2)./sinc((u-uk(k))/2);
end

plot(u,Bc,'g')

Bk = conj(Bp).*exp(-j*pi*uk*(N-1)/2);
b = ifft(Bk);
w = b.*exp(-j*[0:1:N-1]*pi*((N-1)/2-delta)*2/N);
w = w.';
Bw = w'*v;
plot(u,real(Bw),'r')
xlabel('u=\psi/\pi')
ylabel('Beampattern')
title(['Problem 3.3.2, N = ' int2str(N) ', \Delta = ' num2str(delta)])
hold off

subplot(2,1,2)
hold on
plot(u,20*log10(abs(max(Bd,1e-3))))
plot(uk,20*log10(abs(max(Bp,1e-3))),'o')
plot(u,20*log10(abs(Bc)),'g')
plot(u,20*log10(abs(Bw)),'r')
axis([-1 1 -60 0])
hold off
xlabel('u=\psi/\pi')
ylabel('Beampattern (dB)')

set(gcf,'Paperposition',[0.25 1 8 9])
