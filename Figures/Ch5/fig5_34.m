%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.34
% Wavenumber spectra for complex AR(2) process with additive white noise
% sigmau^2/sigmaw^2 = 0 dB, z1=0.5, z2 = 0.5exp(j0.7pi)
% K. Bell 7/24/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

sgw_u = 10^(0/10);
psi = [-1:0.0005:1]*pi;
nf = size(psi,2);

%(a)
z1 = 0.5;
z2 = 0.5*exp(j*0.7*pi);

a1 = -(z1+z2);
a2 = z1*z2;

% find roots of 
% 1+sgw_u*abs(1+a1*z^-1 + a2*z^-2)^2= 1+sgw_u(conj(p2)*z^2 +conj(p1)*z +p0 +p1*z^-1 +p2*z^-2)
% = sgb_u*abs((1-zb1*z^-1)*(1-zb2*z^-1))^2
% = sgb_u*abs((1+b1*z^-1 + b2*z^-2)^2
% Then b1 = -(zb1+zb2)
% b2 = zb1*zb2
% sgb_u = 1+sgw_u*abs(1+a1+a2)^2 / abs(1+b1+b2)^2

p2 = sgw_u*a2;
p1 = sgw_u*(a1+conj(a1)*a2);
p0 = sgw_u*(1+abs(a1)^2 +abs(a2)^2);

p = [conj(p2) conj(p1) 1+p0 p1 p2];
z = roots(p);
[y,I]=sort(abs(z));
zb1 = z(I(1));
zb2 = z(I(2));
b1 = -(zb1+zb2);
b2 = zb1*zb2;
sgb_u = (1+sgw_u*(abs(1+a1+a2)^2))/(abs(1+b1+b2)^2);

P = (ones(1,nf)./((abs(ones(1,nf)+a1*exp(-j*psi)+a2*exp(-j*2*psi)) ).^2)+ sgw_u*ones(1,nf));     % actual
Pb = sgb_u*((abs(ones(1,nf)+b1*exp(-j*psi)+b2*exp(-j*2*psi)) ).^2)./((abs(ones(1,nf)+a1*exp(-j*psi)+a2*exp(-j*2*psi)) ).^2);  %ARMA(2,2)
plot(psi/pi,10*log10(P))
xlabel('\psi/\pi','Fontsize',14)
ylabel(['{\itP}_{\itxx}(\psi)/\sigma^2_{\itu} (dB)'],'Fontsize',14)
grid on
hold on
plot(psi/pi,10*log10(Pb),'r')

L2 = 6;
c = [a1-b1 -b1*(a1-b1)+a2-b2 zeros(1,L2-2)];
zv = exp(-j*[0:1:L2].'*psi);
Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
plot(psi/pi,10*log10(Pc),'--r')
for m=3:L2
    c(m) = -c(m-1)*b1 - c(m-2)*b2; 
    Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
end
plot(psi/pi,10*log10(Pc),'-.r')
hold off
h=legend('Actual','AR(2)',['AR(' int2str(L2) ')'],2);
set(h,'Fontsize',12)
