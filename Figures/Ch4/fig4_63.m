%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.63
% Response of cylindrical array versus theta for phi = 0 and 0.045*pi
% Xin Zhang 9/21/99
% Updated by K. Bell 10/12/00
% Updated by Lillian Xu 12/06/2000, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters:
%    N=25

clear all
close all

Rl=5/pi;
N=25;
Ms=21;
M=(Ms-1)/2;
phi_set=[0 0.045]*pi;

Mm=(-M:1:M)';
theta=0:0.0025*pi:pi;

nbar = 6;
d = Rl*2/N;
R = 10^(20/20);
x0 = cosh(acosh(R)/(Ms-1));
cc=cos((2*[1:1:floor((Ms-1)/2)]-1)*0.5*pi/(Ms-1)).';
udc = acos(cc/x0)/(pi*d);

if rem(N,2)==0
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

Nz = 11;
Dz = (-(Nz-1)/2:1:(Nz-1)/2).';

R = 10^(30/20);
x0 = cosh(acosh(R)/(Nz-1));
wt = poly(exp(j*2*acos(cos((2*[1:1:Nz-1]-1)*pi/(2*(Nz-1)))/x0))).';
wz = real(wt/sum(wt));
vz = exp(j*Dz*pi*cos(theta));
Bd = wz'*vz;

Bt = 1;
for p = 1:length(phi_set)
   phi = phi_set(p);
   for l = 1:length(theta)
   
for n=1:Ms
   m=n-1-M;
   if m<0
      Vj(n,1)=besselj(abs(m),2*pi*Rl*sin(theta(l)))*(-1)^m;
   else
      Vj(n,1)=besselj(m,2*pi*Rl*sin(theta(l)));
   end
end
J=diag(Vj);
V=exp(j*Mm*phi);
B1(l)=Wpm'*J*V;
end

B = B1.*Bd;
if p == 1
   Bt = max(abs(B));
end
B = abs(B)/Bt;
Bm=20*log10(abs(B));

figure
plot(theta/pi,Bm)
grid
axis([0 1 -60 0])
xlabel('\theta/\pi','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
%if p == 1
%title('Response of Cylindrical Array, \phi = 0')
%else
%title('Response of Cylindrical Array, \phi = 0.045\pi')
%end

end
