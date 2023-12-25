%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.2
% Raised cosine weighting function
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
p=[0;0.17;0.31;1];
W_rcos = p*ones(1,M)+(ones(4,1)-p)*cos(pi*D*2/M);

% Beampatterns
u = [0:0.001:1];
c=p+(ones(4,1)-p)*2/pi;
A = exp(-j*2*pi*D.'*u);
G_rcos = W_rcos*A;
for m=1:4
 G_rcos(m,:) = G_rcos(m,:)/(max(abs(G_rcos(m,:))));
end

figure
clf
plot(u,20*log10(abs(G_rcos(3,:))),'-.')
hold on
plot(u,20*log10(abs(G_rcos(2,:))),'--')
plot(u,20*log10(abs(G_rcos(1,:))),'-')
hold off
axis([0 1 -80 0])
h=legend('{\it p}=0.31','{\it p}=0.17','{\it p}=0');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
