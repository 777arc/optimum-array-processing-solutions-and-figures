%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3.2.6(a)
% K. Bell 10/14/98
% Updated 9/19/00 by K. Bell
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=15;
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
r=zeros((N-1)/2,3);

figure(1)

for n=1:3
   subplot(2,2,n)
   z = roots(wn(:,n));
   abs(z)
   a = sort(angle(z));
   r(:,n)=a((N+1)/2:N-1)/pi;
   plot(real(z),imag(z),'o')
   hold on
   plot(x,y)
   hold off
   title(['DPSS, \psi_o = ' num2str(psi0(n)) '\pi, N = ' int2str(N)])
   axis('square')
   axis([-1.2 1.2 -1.2 1.2])
   grid on
end


set(gcf,'Paperposition',[0.25 1 8 9])

r
