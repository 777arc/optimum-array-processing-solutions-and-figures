%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.15 & 2.16
% Beampattern of uniformly weighted linear array
% Xiaomin Lu
% Updated 1/5/99
% Last updated  by K. Bell 7/22/01, 10/4/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 11;				% 11 elements
% Note, Eq. (2.92) is in fact the standard Dirichlet function
% beam = diric(psi, N);
% if you don't have the signal processing toolbox, this does the same thing
psi = (-5:1/400:5)*pi;
beam = zeros(1,length(psi(:)));
y=sin(.5*psi);
i=find(abs(y)>1e-12);            % set where x is not divisible by 2 pi
j=1:length(psi(:));
j(i)=[];                         % complement set
beam(i)=sin((N/2)*psi(i))./(N*y(i));
beam(j)=sign(cos(psi(j)*((N+1)/2)));

plot(psi/pi, beam)
grid
xlabel('{\it \psi}/\pi','Fontsize',14)
ylabel('Frequency wavenumber response function','Fontsize',14)
axis([-5 5 -0.4 1])

figure
beam = 20*log10(abs(beam));
plot(psi/pi,beam);
xlabel('{\it\psi/}\pi','Fontsize',14)
ylabel('Frequency wavenumber response function (dB)','Fontsize',14)
axis([-5 5 -25 0])
grid

