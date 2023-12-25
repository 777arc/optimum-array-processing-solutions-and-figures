%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.23
% Wavenumber spectra for complex AR(1) process
% |a(1)| = 0.5, 0.7, and 0.9; phi_a = 0
% K. Bell 10/13/00
% updated by K. Bell 7/24/01, 10/23/01, 11/15/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

a1= -[0.5 0.7 0.9];
phi=0;
psi = [-1:0.01:1]*pi;
nf = size(psi,2);

n1 = size(a1,2);
P=zeros(n1,nf);

for n=1:n1
   z1 = -a1(n)*exp(-j*pi*phi);
   P(n,:) = ones(1,nf)./((abs(ones(1,nf)-z1*exp(-j*psi)) ).^2);
end

plot(psi/pi,10*log10(P(1,:)))
hold on
plot(psi/pi,10*log10(P(2,:)),'--')
plot(psi/pi,10*log10(P(3,:)),'-.')
h=legend(['|{\itz}_1|=',num2str(a1(1))],['|{\itz}_1|=',num2str(a1(2))],...
   ['|{\itz}_1|=',num2str(a1(3))]);
set(h,'Fontsize',12)
xlabel('\psi/\pi','Fontsize',14)
ylabel(['{\itP}_{\itxx}(\psi)/\sigma^2_{\itu} (dB)'],'Fontsize',14)
grid on