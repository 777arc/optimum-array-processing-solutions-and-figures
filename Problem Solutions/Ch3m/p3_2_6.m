%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3.2.6
% K. Bell 10/14/98
% Updated 9/19/00 by K. Bell
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=13;
d=0.5;
D = d*[-(N-1)/2:1:(N-1)/2].';
u = [-1:0.01:1];
v = exp(j*2*pi*D*u);
w=zeros(N,3);
wn = zeros(N,3);

psi0 = [0.15 0.25 0.5];

for k=1:3
  r = 2*pi*psi0(k)*sinc([0:1:N-1]*psi0(k));
  A = toeplitz(r,r);
  [vv,l]=eig(A);
  [l,ind]=sort(diag(l));
  vv=vv(:,ind);
  ww=vv(:,N);
  wn(:,k)=ww/sum(ww);
end

t=[0:0.01:1]*2*pi;
x=cos(t);
y=sin(t);
r=zeros((N-1)/2,5);
zz=zeros((N-1)/2,5);

figure(1)
subplot(3,2,1)
z = roots(wn(:,1));
a = sort(angle(z));
r(:,1)=a((N+1)/2:N-1)/pi;
plot(real(z),imag(z),'o')
hold on
plot(x,y,'--')
hold off
title(['DPSS, \psi_o = 0.15\pi, N = ' int2str(N)])
axis('square')
axis([-1.2 1.2 -1.2 1.2])
grid on
subplot(3,2,3)
z = roots(wn(:,2));
a= sort(angle(z));
r(:,2)=a((N+1)/2:N-1)/pi;
plot(real(z),imag(z),'o')
hold on
plot(x,y,'--')
hold off
grid on
title(['DPSS, \psi_o = 0.25\pi, N = ' int2str(N)])
axis('square')
axis([-1.2 1.2 -1.2 1.2])

subplot(3,2,5)
z = roots(wn(:,3));
a = sort(angle(z));
r(:,3)=a((N+1)/2:N-1)/pi;
plot(real(z),imag(z),'o')
hold on
plot(x,y,'--')
hold off
grid on
title(['DPSS, \psi_o = 0.5\pi, N = ' int2str(N)])
axis('square')
axis([-1.2 1.2 -1.2 1.2])

% Kaiser
w=zeros(N,2);
wn = zeros(N,2);

beta = [2 8];

for k=1:2
  ww = besseli(0,beta(k)*sqrt(1-(2*(D/d)/(N-1)).^2));
  wn(:,k)=ww/sum(ww);
end


figure(1)
subplot(3,2,2)
z = roots(wn(:,1));
a = sort(angle(z));
r(:,4)=a((N+1)/2:N-1)/pi;
plot(real(z),imag(z),'o')
hold on
plot(x,y,'--')
hold off
title(['Kaiser, \beta = 2, N = ' int2str(N)])
axis('square')
axis([-1.2 1.2 -1.2 1.2])
grid on
subplot(3,2,4)
z = roots(wn(:,2));
a = sort(angle(z));
r(:,5)=a((N-1)/2:N-2)/pi;
plot(real(z),imag(z),'o')
hold on
plot(x,y,'--')
hold off
grid on
title(['Kaiser, \beta = 8, N = ' int2str(N)])
axis('square')
axis([-1.2 1.2 -1.2 1.2])


set(gcf,'Paperposition',[0.25 1 8 9])

r
