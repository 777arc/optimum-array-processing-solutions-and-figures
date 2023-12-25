%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.12
% Zero plots for Hamming, Blackman-Harris weightings
% Lillian Xu 3/18/98
% K. Bell 7/22/01,9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=11; 
us=0;  

u=-1:0.001:1;
n=conj(-(N-1)/2:(N-1)/2)';
vs=exp(i*n.*pi.*us); 

%------------------------Weighting
w1=0.54+0.46*cos(2*pi*n/N);     		% Hamming
w2=0.42+0.5*cos(2*pi*n/N)+0.08*cos(4*pi*n/N);	% Blackman-Haris

z1=roots(w1);
z2=roots(w2);

% ---------------- Zero Plot.
figure

subplot(1,2,1);
plot(real(z1),imag(z1),'o');
grid on;
axis([-1.2 1.2 -1.2 1.2])
axis('square')
title('Hamming','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
hold off

subplot(1,2,2);
plot(real(z2),imag(z2),'o');
grid on;
axis([-1.2 1.2 -1.2 1.2])
axis('square')  
title('Blackman-Harris','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
hold on
theta=2*pi*[0:0.01:1];
plot(cos(theta),sin(theta),'--');
text(-1.2,-1.7,(['Another zero is at (',num2str(real(z2(1))),', ',num2str(imag(z2(1))),')']),'Fontsize',12)
hold off
