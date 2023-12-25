%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.32
% Wavenumber spectra for complex AR(1) process with additive white noise
% sigmau^2/sigmaw^2 = -10 dB, |a1|=0.5
% K. Bell 7/24/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

sgw_u = 10^(10/10);
psi = [-1:0.0005:1]*pi;
nf = size(psi,2);

%(a)
a1 = -0.5;
x = 1+sgw_u+sgw_u*abs(a1)^2;
y = 4*(sgw_u*abs(a1))^2;
sgb_u = 0.5*(x+sqrt(x^2-y));

b1 = (sgw_u/sgb_u)*a1;
P = (ones(1,nf)./((abs(ones(1,nf)+a1*exp(-j*psi)) ).^2)+ sgw_u*ones(1,nf));     % actual
Pb = sgb_u*((abs(ones(1,nf)+b1*exp(-j*psi)) ).^2)./((abs(ones(1,nf)+a1*exp(-j*psi)) ).^2);  %ARMA 1
plot(psi/pi,10*log10(P))
xlabel('\psi/\pi','Fontsize',14)
ylabel(['{\itP}_{\itxx}(\psi)/\sigma^2_{\itu} (dB)'],'Fontsize',14)
grid on
hold on
plot(psi/pi,10*log10(Pb),'r')
L1 = 2;
L2 = 6;
c = [a1-b1 zeros(1,L2-1)];
zv = exp(-j*[0:1:L2].'*psi);
Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
plot(psi/pi,10*log10(Pc),'--r')
for m=2:L1
    c(m) = -c(m-1)*b1; 
    Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
end
plot(psi/pi,10*log10(Pc),':r')
for m=L1+1:L2
    c(m) = -c(m-1)*b1; 
    Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
end
plot(psi/pi,10*log10(Pc),'-.r')
hold off
h=legend('Actual','AR(1)',['AR(' int2str(L1) ')'],['AR(' int2str(L2) ')']);
set(h,'Fontsize',12)


%(b)
a1 = -0.9;
x = 1+sgw_u+sgw_u*abs(a1)^2;
y = 4*(sgw_u*abs(a1))^2;
sgb_u = 0.5*(x+sqrt(x^2-y));

b1 = (sgw_u/sgb_u)*a1;
P = (ones(1,nf)./((abs(ones(1,nf)+a1*exp(-j*psi)) ).^2)+ sgw_u*ones(1,nf));     % actual
Pb = sgb_u*((abs(ones(1,nf)+b1*exp(-j*psi)) ).^2)./((abs(ones(1,nf)+a1*exp(-j*psi)) ).^2);  %ARMA 1
figure
plot(psi/pi,10*log10(P))
xlabel('\psi/\pi','Fontsize',14)
ylabel(['{\itP}_{\itxx}(\psi)/\sigma^2_{\itu} (dB)'],'Fontsize',14)
grid on
hold on
plot(psi/pi,10*log10(Pb),'r')
L1 = 4;
L2 = 16;
c = [a1-b1 zeros(1,L2-1)];
zv = exp(-j*[0:1:L2].'*psi);
Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
plot(psi/pi,10*log10(Pc),'--r')
for m=2:L1
    c(m) = -c(m-1)*b1; 
    Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
end
plot(psi/pi,10*log10(Pc),':r')
for m=L1+1:L2
    c(m) = -c(m-1)*b1; 
    Pc = sgb_u*ones(1,nf)./(abs([1 c]*zv).^2);
end
plot(psi/pi,10*log10(Pc),'-.r')
hold off
h=legend('Actual','AR(1)',['AR(' int2str(L1) ')'],['AR(' int2str(L2) ')']);
set(h,'Fontsize',12)

