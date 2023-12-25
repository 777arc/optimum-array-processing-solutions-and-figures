%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.48
% K. Bell 9/29/00, 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 10;                         % number of elements            
d = 0.5;                        % spacing wrt lambda
D = [-(N-1)/2:1:(N-1)/2]'*d;     % positions 

% steering direction
uT = 0 ;                         
vT = exp(j*2*pi*D*uT);

uu = [-0.2:0.01:0.2];
vv = exp(j*2*pi*D*uu);

figure(1)
plot([-1 uu 1],[0 abs(vT'*vv)/N 0])
axis([-1 1 0 1])
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)
