%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.1
% Cosine weighting: Array beam patterns
% Kristine Bell
% Modified by Xin Zhang 3/17/99
% Lillian 04/16/2001, K. Bell 7/20/01, K. Bell 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

M=11;
d=0.5;                        % sensor spacing wrt wavelength
D = [-(M-1)/2:1:(M-1)/2]*d;   % sensor positions in wavelengths
% weights, normalized so that w(0)=1
W_unf = ones(1,M);
W_cos = cos(pi*D*2/M);

% Beampatterns 
u = [0:0.001:1];              
A = exp(-j*2*pi*D.'*u);
G_unf = W_unf*A;
G_unf = G_unf/(max(abs(G_unf)));
G_cos = W_cos*A;
G_cos = G_cos/(max(abs(G_cos)));
figure
% array only
clf
plot(u,20*log10(abs(G_unf)),'--')
hold on
plot(u,20*log10(abs(G_cos)),'-')
hold off
axis([0 1 -80 0])
%title([int2str(M) ' element array'])
h=legend('Uniform','Cosine');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
