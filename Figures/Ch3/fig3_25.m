%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.25
% Chebychev polynomial of the seventh degree
% Lillian Xu
% Xin Zhang 3/23/99, Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=8;
x=-1.2:0.001:1.2;

sidelobe=26;   	% *dB below the main lobe maximum
R=10^(sidelobe/20);
x0_set=cosh(1/(N-1)*acosh(R));

%x=x0*cos(pi*u/2);

T7=64*x.^7-112*x.^5+56*x.^3-7*x;

figure
plot(x,T7);
hold on
plot([-1.2 1.2],[1 1],'--');
plot([-1.2 1.2],[-1 -1],'--');
plot([x0_set x0_set],[-24 24],'--');
xlabel('\it x','Fontsize',14)
ylabel('{\it T}_7{\it(x)}','Fontsize',14)
%title(['Dolph-Chebychev Polynomial'])
grid on
axis([-1.5 1.5 -30 30])
text(x0_set-0.02,-26,'{\it x}_{0}','Fontsize',12)