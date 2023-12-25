%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.24
% Wavenumber spectra for complex AR(1) process
% |a(1)| = 0.7; phi_a = 0, 0.6, and 1.0
% K. Bell 10/13/00
% updated by K. Bell 7/24/01,10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

a1= -0.7;
phi=[0.0 0.6 1];
psi = [-1:0.01:1]*pi;
nf = size(psi,2);


n1 = size(phi,2);
P=zeros(n1,nf);

for n=1:n1
   z1 = -a1*exp(j*pi*phi(n));
   P(n,:) = ones(1,nf)./((abs(ones(1,nf)-z1*exp(-j*psi)) ).^2);
end

plot(psi/pi,10*log10(P(1,:)))
hold on
plot(psi/pi,10*log10(P(2,:)),'--')
plot(psi/pi,10*log10(P(3,:)),'-.')
h=legend(['\phi=',num2str(phi(1))],['\phi=',num2str(phi(2))],...
   ['\phi=',num2str(phi(3))]);
set(h,'Fontsize',12)
xlabel('\psi/\pi','Fontsize',14)
ylabel(['{\itP}_{\itxx}(\psi)/\sigma^2_{\itu} (dB)'],'Fontsize',14)

grid on
