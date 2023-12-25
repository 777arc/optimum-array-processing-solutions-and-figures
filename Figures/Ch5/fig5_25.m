%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.25
% Wavenumber spectra for complex AR(2) process
% z_1 = 0.9, z_2 = 0.7exp(j0.8pi)
% K. Bell 10/13/00
% updated by K. Bell 7/24/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

psi = [-1:0.001:1]*pi;
nf = size(psi,2);

z1 = 0.9*exp(j*pi*0);
z2 = 0.7*exp(j*pi*0.8);

P = ones(1,nf)./(((abs(ones(1,nf)-z1*exp(-j*psi))).^2).*((abs(ones(1,nf)-z2*exp(-j*psi))).^2));


plot(psi/pi,10*log10(P))

xlabel('\psi/\pi','Fontsize',14)
ylabel(['{\itP}_{\itxx}(\psi)/\sigma^2_{\itu} (dB)'],'Fontsize',14)
grid on