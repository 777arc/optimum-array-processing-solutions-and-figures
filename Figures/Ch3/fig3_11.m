%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.11
% Zero plots for cosine, cosine squared weightings
% Lillian Xu
% Modified by Xin Zhang 3/18/98
% K. Bell 7/20/01, 9/30/01, 10/17/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=11; 
n=conj(-(N-1)/2:(N-1)/2)';

w2=cos(pi*n/N);
w3=w2.*w2;
z2=roots(w2);
z3=roots(w3);

% ---------------- Zero Plot.
figure
%%%%%%%%%%%%	plot 1
subplot(1,2,1);
plot(real(z2),imag(z2),'o');
axis([-1.2 1.2 -1.2 1.2])
axis('square')
grid on;
title('Cosine','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
plot([-0.9 -0.65],[0.03 0.17])
text(-0.6,0.2,'2 zeros','Fontsize',12)
hold off

%%%%%%%%%%%	plot 2
subplot(1,2,2);
plot(real(z3),imag(z3),'o');
grid on;
axis([-1.2 1.2 -1.2 1.2])
axis('square')
title('Cosine-squared','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
text(-1.2,-1.7,(['Another zero is at (',num2str(real(z3(1))),', ',num2str(imag(z3(1))),')']),'Fontsize',12)
%text(1.6,-0.5,(['another zero at :']));
%text(1.7,-1,(['( ',num2str(real(z3(1))),' , ',num2str(imag(z3(1))),' )']))
hold off
