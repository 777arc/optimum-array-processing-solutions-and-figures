%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.43
% Beampattern of uniformly weighted circular array
% Xiaomin Lu 11-2-98	
% Updated by K. Bell 10/2/00		
% Updated by Lillian Xu 12/06/2000, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

u = 0:1/200:20;
beam = besselj(1,pi*u)./(pi*u);
beam = 20*log10(abs(beam)/max(abs(beam)));
plot(u,beam)
grid
axis([0 20 -60 0])
xlabel('\it u_R','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
%title('Beampattern of uniformly weighted circular array')
