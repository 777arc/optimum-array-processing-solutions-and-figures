%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.33
% Comparison of J0(psi_R) and sinc(psi_R)
% Xin Zhang  9/20/99	
% Updated by K. Bell 9/29/00, 9/30/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

psi_R = 0:0.1:10;
J0 = besselj(0,psi_R);
SC = sinc(psi_R/pi);

figure
plot(psi_R,J0,'-',psi_R,SC,'--')
grid on
axis([0 6.25 -0.5 1])
h=legend('{\it J_0(\psi_R)}','sinc{\it(\psi_R)}');
set(h,'Fontsize',12)
xlabel('\it\psi_R','Fontsize',14)
