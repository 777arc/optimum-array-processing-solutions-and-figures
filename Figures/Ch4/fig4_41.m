%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.41
% Response of cylindrical array versus phi for theta = 0.5137*pi
% Xin Zhang 9/21/99
% Updated by K. Bell 9/30/00
% Updated by Lillian Xu 12/06/2000, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Description:
%    Redrawing of Figure 4.41,
% Parameters:
%    N=25

clear all
close all

Rl=5/pi;
N=25;
Ms=21;
M=(Ms-1)/2;
theta=0.5137*pi;

Mm=(-M:1:M)';
phi=(-1)*pi:0.002*pi:pi;

SL = -20;
nbar = 6;
d = Rl*2/N;
R = 10^(-SL/20);
x0 = cosh(acosh(R)/(Ms-1));
cc=cos((2*[1:1:floor((Ms-1)/2)]-1)*0.5*pi/(Ms-1)).';
udc = acos(cc/x0)/(pi*d);

if rem(Ms,2)==0
   uu = [1:1:(Ms/2)-1].'/(Ms*d)
   u_nbar = udc(nbar);       % u_nbar Chebychev
   sigma = nbar/(Ms*d*u_nbar);
   udcmod = udc*sigma;           % modified Chebychev
   uz = [udcmod(1:nbar-1);uu(nbar:(Ms/2)-1)];
   ut = [uz;-uz;-1]
else
   uu = [1:1:(Ms-1)/2].'/(Ms*d);
   u_nbar = udc(nbar) ;      % u_nbar Chebychev
   sigma = nbar/(Ms*d*u_nbar);
   udcmod = udc*sigma;           % modified Chebychev
   uz = [udcmod(1:nbar-1);uu(nbar:(Ms-1)/2)];
   ut = [uz;-uz];
end
ut;
w = poly(exp(j*2*pi*d*ut)).';
Wh=w/sum(w);

for n=1:Ms
   m=n-1-M;
   if m<0
      Wpm(n,1)=Wh(n,1)/besselj(abs(m),2*pi*Rl)/(-1)^m;
   else
      Wpm(n,1)=Wh(n,1)/besselj(m,2*pi*Rl);
   end
end
Wpm = Wpm./sum(Wpm);
for n=1:Ms
   m=n-1-M;
   if m<0
      Vj(n,1)=besselj(abs(m),2*pi*Rl*sin(theta))*(-1)^m;
   else
      Vj(n,1)=besselj(m,2*pi*Rl*sin(theta));
   end
end
J=diag(Vj);
V=exp(j*Mm*phi);
B1=Wpm'*J*V;
B1 = abs(B1)/max(abs(B1));

Ll = 10;
Bd = sinc(Ll*cos(theta));
B = B1.*Bd;
Bm=20*log10(abs(B));

figure
plot(phi/pi,Bm)
grid
axis([-1 1 -60 0])
xlabel('\phi/\pi','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
%title('Response of Cylindrical Array, \theta = 0.5137\pi')