%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.4.13
% K. Bell 10/19/99
% Last Updated 9/02/02 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=41;
nbar = [6 10 15];
SL = -30;
m=N-1;
d=0.5;
D = d*[-(N-1)/2:1:(N-1)/2].';
u = [-1:0.001:1];
v = exp(j*2*pi*D*u);
nu = length(u);

wn = zeros(N,3);

R = 10^(-SL/20);
x0 = cosh(acosh(R)/m);
cc=cos((2*[1:1:floor((N-1)/2)]-1)*0.5*pi/(N-1)).';   % Dolph-chebychev roots
udc = acos(cc/x0)/(pi*d);

for qq=1:3
   if rem(N,2)==0                   % N even
      uu = [1:1:(N/2)-1].'/(N*d)    % uniform roots
      u_nbar = udc(nbar(qq));       % u_nbar Chebychev
      sigma = nbar(qq)/(N*d*u_nbar);
      udcmod = udc*sigma;           % modified Chebychev
      uz = [udcmod(1:nbar(qq)-1);uu(nbar(qq):(N/2)-1)];
      ut = [uz;-uz;-1];
   else                             % Nodd
      uu = [1:1:(N-1)/2].'/(N*d);   % uniform roots
      u_nbar = udc(nbar(qq)) ;      % u_nbar Chebychev
      sigma = nbar(qq)/(N*d*u_nbar);
      udcmod = udc*sigma;           % modified Chebychev
      uz = [udcmod(1:nbar(qq)-1);uu(nbar(qq):(N-1)/2)];
      ut = [uz;-uz];
   end
   
   w = poly(exp(j*2*pi*d*ut)).';
   wn(:,qq)=w/sum(w);
   
   Bw = wn(:,qq)'*v;
   figure(1)
   subplot(3,1,qq)
   plot(u,10*log10(abs(Bw).^2),'r')
   grid
   xlabel('u=\psi/\pi')
   ylabel('Beampattern (dB)')
   set(gca,'YTick',[-60:10:0])   
   if qq ==1
      title(['Problem 3.4.13, Villeneuve N = ' int2str(N) ', SL = ' num2str(SL) ', nbar = ' int2str(nbar(qq))])
   else
      title(['nbar = ' int2str(nbar(qq))])
   end
   axis([-1 1 -60 0])
   
end
set(gcf,'Paperposition',[0.25 1 8 9])

wn(:,1)
