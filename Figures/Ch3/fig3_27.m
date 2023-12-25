%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.27
% Zero plots of the Dolph-Chebychev weightings
%   (a) -20 dB
%   (b) -30 dB
%   (c) -40 dB
% Lillian Xiaolan Xu 3/23/99
% K. Bell 7/22/01, 9/30/01, 10/17/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

N=8;
sidelobe=[20 30 40];   	% *dB below the main lobe maximum
R=10.^(sidelobe/20);
x0_set=cosh(1/(N-1)*acosh(R));

for num=1:size(x0_set,2)
    x0=x0_set(num);
    a3=x0^7;
    a2=(a3-x0^5)*112/16;
    a1=((x0^3-a3)*56+20*a2)/4;
    a0=-7*x0+3*a1-5*a2+7*a3;
    
    wdq=[a3 a2 a1 a0 a0 a1 a2 a3];        %w'
    sumw=sum(wdq);
    wdq=wdq./sumw;
    wdq=wdq';
    W(:,num)=wdq;
    
    %x=x0*cos(pi*u/2);
    %T7=64*x.^7-112*x.^5+56*x.^3-7*x;
end;

z1=roots(W(:,1));
z2=roots(W(:,2));
z3=roots(W(:,3));

%-----------------------------------------------------------%
% ---------------- Zero Plot.
figure
%%%%%%%%%%%%      plot 1

subplot(2,2,1);
plot(real(z1),imag(z1),'o');
axis([-1.2 1.2 -1.2 1.2])
axis('square')
grid on;
title('-20 dB sidelobes','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
hold off

%%%%%%%%%%%%	plot 2

subplot(2,2,2);
plot(real(z2),imag(z2),'o');
axis([-1.2 1.2 -1.2 1.2])
axis('square')
grid on;
title('-30 dB sidelobes','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
hold off

%%%%%%%%%%%	plot 3

subplot(2,1,2);
plot(real(z3),imag(z3),'o');
grid on;
axis([-1.2 1.2 -1.2 1.2])
axis('square')
title('-40 dB sidelobes','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
hold off
