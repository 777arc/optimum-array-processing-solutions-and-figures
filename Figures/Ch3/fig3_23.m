%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.23
% Chebychev polynomials:
%      (a) n=2
%      (b) n=3
%      (c) n=4
%
% Lillian Xu 04/16/2001
% updated by K. Bell 7/22/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

x=-1.2:0.001:1.2;


T2=2*x.^2-1;
T3=4*x.^3-3*x;
T4=8*x.^4-8*x.^2+1;

figure

subplot(1,3,1)
plot(x,T2);
hold on
plot([-1.2 1.2],[1 1],'--');
plot([-1.2 1.2],[-1 -1],'--');
xlabel('\it x','Fontsize',14)
ylabel('{\it T}_2{\it(x)}','Fontsize',14)

grid on


subplot(1,3,2)
plot(x,T3);
hold on
plot([-1.2 1.2],[1 1],'--');
plot([-1.2 1.2],[-1 -1],'--');
xlabel('\it x','Fontsize',14)
ylabel('{\it T}_3{\it(x)}','Fontsize',14)
%title(['Dolph-Chebychev Polynomial'])
grid on



subplot(1,3,3)
plot(x,T4);
hold on
plot([-1.2 1.2],[1 1],'--');
plot([-1.2 1.2],[-1 -1],'--');
xlabel('\it x','Fontsize',14)
ylabel('{\it T}_4{\it(x)}','Fontsize',14)

grid on



