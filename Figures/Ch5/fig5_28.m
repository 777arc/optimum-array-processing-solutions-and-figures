%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.28
% Wavenumber spectra for complex AR(3) process
% z_1 = 0.9exp(-j0.7pi), z_2 = 0.9exp(-j0.5pi), z_3 = 0.7exp(-0.4pi)
% K. Bell 10/13/00
% updated by K. Bell 7/24/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

psi = [-1:0.0005:1]*pi;
nf = size(psi,2);

z1 = 0.9*exp(-j*pi*0.7);
z2 = 0.9*exp(-j*pi*0.5);
z3 = 0.7*exp(j*pi*0.4);

P = ones(1,nf)./(((abs(ones(1,nf)-z1*exp(-j*psi)) ).^2).*...
   ((abs(ones(1,nf)-z2*exp(-j*psi))).^2).*((abs(ones(1,nf)-z3*exp(-j*psi))).^2));

plot(psi/pi,10*log10(P))
xlabel('\psi/\pi','Fontsize',14)
ylabel(['{\itP}_{\itxx}(\psi)/\sigma^2_{\itu} (dB)'],'Fontsize',14)
grid on