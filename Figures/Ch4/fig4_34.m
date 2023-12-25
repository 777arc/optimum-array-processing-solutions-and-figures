%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.34
% J0(psi_R) through J7(psi_R) vs. psi_R; psi_R = 2*pi*R_lambda*sin(theta): R_lambda = 1
% Xin Zhang 9/20/99	
% Updated by K. Bell 9/29/00, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

psi_R = 0:0.1:10;
J0 = besselj(0,psi_R);
J1 = besselj(1,psi_R);
J2 = besselj(2,psi_R);
J3 = besselj(3,psi_R);
J4 = besselj(4,psi_R);
J5 = besselj(5,psi_R);
J6 = besselj(6,psi_R);
J7 = besselj(7,psi_R);

figure
plot(psi_R,J0,'-')
hold on;
plot(psi_R,J1,'--')
plot(psi_R,J2,'-.')
plot(psi_R,J3,':')
plot(psi_R,J4,'-')
plot(psi_R,J5,'--')
plot(psi_R,J6,'-.')
plot(psi_R,J7,':')
hold off;
grid on
axis([0 6.25 -0.5 1.5])
h=legend('\itJ_0(\psi_R)','J_1(\psi_R)','J_2(\psi_R)','J_3(\psi_R)', ...
   'J_4(\psi_R)','J_5(\psi_R)','J_6(\psi_R)','J_7(\psi_R)');
set(h,'Fontsize',12)
xlabel('\it\psi_R','Fontsize',14)
