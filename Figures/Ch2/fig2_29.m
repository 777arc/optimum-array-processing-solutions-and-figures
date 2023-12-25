%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.29
% Beam pattern of uniformly weighted linear aperture
% Kristine Bell 1/20/99
% Updated by K. Bell 6/25/01, 10/4/01
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

arg= [-5:0.01:5];       % u*L/lambda
na=length(arg);
G0   = sinc(arg);
h=plot(arg,G0);
%set(h,'LineWidth',1.5)
hold on
plot([-5 5],[0 0])
set(gca,'Xtick',[-5:1:5])
grid on
hold off
ylabel('Beam pattern')
xlabel('{\ituL}/\lambda')


