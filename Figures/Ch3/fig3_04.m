%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.4
% Beam patterns for Hann, Hamming, and Blackman-Harris
% Kristine Bell
% Modified by Xin Zhang 3/17/99
% Lillian Xu 04/16/2001, K. Bell 7/20/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hann, Hamming, Blackman-Harris
% M-element array, d=lambda/2

clear all
close all

M=11;
d=0.5;                        % sensor spacing wrt wavelength
D = [-(M-1)/2:1:(M-1)/2]*d;   % sensor positions in wavelengths

% weights, normalized so that w(0)=1
W_Hann = 0.5*(ones(1,M)+cos(2*pi*D*2/M));
W_Hamm = 0.54*ones(1,M)+0.46*cos(2*pi*D*2/M);
W_BH   = 0.42*ones(1,M)+0.5*cos(2*pi*D*2/M)+0.08*cos(4*pi*D*2/M);

% Beampatterns

u = [0:0.001:1];
A = exp(-j*2*pi*D.'*u);
G_Hann = W_Hann*A;
G_Hann = G_Hann/(max(abs(G_Hann)));
G_Hamm = W_Hamm*A;
G_Hamm = G_Hamm/(max(abs(G_Hamm)));
G_BH = W_BH*A;
G_BH = G_BH/(max(abs(G_BH)));

figure
clf
plot(u,20*log10(abs(G_Hann)),'--')
hold on
plot(u,20*log10(abs(G_Hamm)),'-')
plot(u,20*log10(abs(G_BH)),'-.')
hold off
axis([0 1 -80 0])
h=legend('Hann','Hamming','Blackman-Harris');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
