%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.17
% Polar plot of B (Theta)
% Lillian Xu 1/5/99
% Updated by K. Bell 6/25/01
% Functions called: polardb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all 
close all

N = 11;                                % Elements in array
d = 0.5;                               % spacing wrt wavelength
beamwidth = 2/(N*d);                   % null-to-null (LL BW is half this)
D=d*[-(N-1)/2:1:(N-1)/2];              % element locations
ang = pi*[-1:0.001:1];
u   = cos(ang); 
n2=size(u,2);
AS  = exp(-j*2*pi*cos(90/180*pi)*D');    % BP points to 90 (broadside)
Au  = exp(-j*2*pi*D'*u);
B   = real(AS'*Au)/N;
G    = 20*log10((abs(B)));
figure(1)
h=polardb(ang,G,-40);
hold off
%set(gca,'Visible','off')
