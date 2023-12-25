%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3.1.4
% K. Bell 10/6/99
% Updated 9/13/00 by K. Bell
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

u = [-1:0.005:1];
u2 = [0:0.0001:0.5];

N=9;
d=0.5;
D = d*[-(N-1)/2:1:(N-1)/2].';
u = [-1:0.01:1];
v = exp(j*2*pi*D*u);
w=zeros(N,4);
wn = zeros(N,4);

psi0 = [0.1 0.2 0.3 0.4];

for k=1:4
  r = 2*pi*psi0(k)*sinc([0:1:N-1]*psi0(k));
  A = toeplitz(r,r);
  [vv,l]=eig(A);
  [l,ind]=sort(diag(l));
  vv=vv(:,ind);
  ww=vv(:,N);
  wn(:,k)=ww/sum(ww);
  w(:,k)=ww/ww((N+1)/2);
end

B = wn'*v;

% Fine grid for computing HPBW, etc. numerically
uu = [0:0.0001:0.5];
mu = length(uu);
v = exp(j*2*pi*D*uu);
BB = wn'*v;

HPBW = zeros(1,4);
BWNN = zeros(1,4);
HPBWdeg = zeros(1,4);
BWNNdeg = zeros(1,4);
SLHT = zeros(1,4);

for k=1:4
   [y,i]=find(abs(BB(k,:))<=sqrt(0.5));
   HPBW(k) = 2*uu(min(i));
   HPBWdeg(k) = 2*(90-acos(HPBW(k)/2)*180/pi);
   
   [y,i] = find(abs(BB(k,2:mu-1))<abs(BB(k,1:mu-2))&abs(BB(k,2:mu-1))<abs(BB(k,3:mu)));
   imin = min(i);
   NULL = uu(imin);
   Ndeg = 90-acos(NULL)*180/pi;
   BWNN(k) = 2*NULL;
   BWNNdeg(k) = 2*Ndeg;

   [y,i]=max(abs(BB(k,imin:mu)));
   SLHT(k) = 10*log10(y^2);
end

WGHTstr  = 'weights    =';
   for m = 2:N
      WGHTstr = [WGHTstr;'            '];
   end
HPBWustr = 'HPBW (u)   =';
HPBWdstr = 'HPBW (deg) =';
BWNNustr = 'BWNN (u)   =';
BWNNdstr = 'BWNN (deg) =';
SLHTstr  = 'SL Ht.(dB) =';
for n=1:4
   HPBWustr = [HPBWustr '  ' sprintf('%8.4f',HPBW(n))];
   HPBWdstr = [HPBWdstr '  ' sprintf('%8.4f',HPBWdeg(n))];
   BWNNustr = [BWNNustr '  ' sprintf('%8.4f',BWNN(n))];
   BWNNdstr = [BWNNdstr '  ' sprintf('%8.4f',BWNNdeg(n))];
   SLHTstr = [SLHTstr '  ' sprintf('%8.4f',SLHT(n))];
   wstr = ['  ' sprintf('%8.4f',w(1,n))];
   for m = 2:N
      wstr = [wstr;['  ' sprintf('%8.4f',w(m,n))]];
   end
   WGHTstr = [WGHTstr wstr];
end
line = ['                                                    '];
list = ['psi_o      =    0.1pi     0.2pi     0.3pi     0.4pi ';...
      line;WGHTstr;line;HPBWustr;HPBWdstr;BWNNustr;BWNNdstr;SLHTstr]

figure(1)
h1=plot(u, 10*log10(abs(B(1,:)).^2),'-');
hold on
h2=plot(u, 10*log10(abs(B(2,:)).^2),'--');
h3=plot(u, 10*log10(abs(B(3,:)).^2),'-.');
h4=plot(u, 10*log10(abs(B(4,:)).^2),':');
h = [h1; h2; h3; h4;];
legend(h,'\psi_o = 0.1\pi','\psi_o = 0.2\pi','\psi_o = 0.3\pi','\psi_o = 0.4\pi')

axis([-1 1 -60 0])
ylabel('Beampattern (dB)')
xlabel('u')
grid on

set(gcf,'Paperposition',[0.25 4 8 6])

