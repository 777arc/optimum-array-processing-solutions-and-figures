%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%% Fig 4.49
%%%% Aperture weighting of J1(0.586p) p =0:pi
%%%% Xiaomin Lu  11-2-98	
%%%% Updated by K. Bell 10/2/00, 9/30/01								
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%


clear all
close all

p = 0:pi/200:pi;
J1 = besselj(1,p);
plot(p,J1)
grid
%title('Aperture Distribution for the Generic Difference Pattern')
xlabel('\it p','Fontsize',14)
ylabel('\it J_1(p)','Fontsize',14)
axis([0 3 0 0.7])


