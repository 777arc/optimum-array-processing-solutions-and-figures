%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.3
% Cosine^m weighting: m = 2,3,and 4
% Kristine Bell
% Modified by Xin Zhang 3/17/99
% Lillian Xu 04/16/2001, K. Bell 7/20/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

M=11;
d=0.5;                        % sensor spacing wrt wavelength
D = [-(M-1)/2:1:(M-1)/2]*d;   % sensor positions in wavelengths

% weights, normalized so that w(0)=1
W_cos2 = cos(pi*D*2/M).^2;
W_cos3 = cos(pi*D*2/M).^3;
W_cos4 = cos(pi*D*2/M).^4;

% Beampatterns

u = [0:0.001:1];
A = exp(-j*2*pi*D.'*u);
G_cos2 = W_cos2*A;
G_cos2 = G_cos2/(max(abs(G_cos2)));
G_cos3 = W_cos3*A;
G_cos3 = G_cos3/(max(abs(G_cos3)));
G_cos4 = W_cos4*A;
G_cos4 = G_cos4/(max(abs(G_cos4)));

figure
clf
plot(u,20*log10(abs(G_cos2)),'--')
hold on
plot(u,20*log10(abs(G_cos3)),'-')
plot(u,20*log10(abs(G_cos4)),'-.')
hold off
axis([0 1 -80 0])
h=legend('{\it m}=2','{\it m}=3','{\it m}=4');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
