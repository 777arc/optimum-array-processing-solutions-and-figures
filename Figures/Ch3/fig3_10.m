%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.10
% Zero plots for 11 element array
% Xin Zhang 3/18/98
% updated by K. Bell 7/20/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 11; 
theta=2*pi*[0:0.01:1];

w = 1/N*ones(N,1);
z = roots(w);
figure
plot(real(z),imag(z),'o');
hold on;
plot(cos(theta),sin(theta),'-'); 
plot(1.5*[-1,1],0*[1,1],'-')
plot(0*[1,1],1.5*[-1,1],'-')
hold off;
set(gca,'YTick',[-1.5 -1 -0.5 0 0.5 1 1.5])
axis square
grid on;
xlabel('{\it x}','Fontsize',14)
ylabel('\it y','Fontsize',14)
