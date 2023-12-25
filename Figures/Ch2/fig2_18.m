%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.18
% beam patterns in different spaces
% Lillian Xiaolan Xu
% Last updated 09/07/2000
% updated by K. Bell 7/22/01, 10/4/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2).';
psi = pi*(-3:0.001:3);
w = 1/N*[ones(N,1)];
d = 1/2;
vv = exp(j*2*d*n*psi);
B = (abs(w'*vv));

u1=-1;
u2=3;
y1=0;
y2=1;

subplot(4,1,1)
plot(-psi/d,B);
hold on
axis([u1*pi/d u2*pi/d y1 y2])
set(gca,'XTick',[-pi/d 0 pi/d 2*pi/d 3*pi/d])
set(gca,'YTick',[0 0.5 1])
set(gca,'XTickLabelMode','manual')
%   set(gca,'XTickLabel',{'-pi/d';'0';'pi/d';'2pi/d';'3pi/d'})
c=0.6;
b=0.15;
text(-pi/d-c,y1-b,'-\pi{\it/d}')
text(0-0.2,y1-b,'0')
text(pi/d-0.4,y1-b,'\pi{\it/d}')
text(2*pi/d-c,y1-b,'2\pi{\it/d}')
text(3*pi/d-c,y1-b,'3\pi{\it/d}')

grid on
legend('{\it k_{z}}-space')
point=0.8;
point_up=0.84;
point_down=0.76;
plot([-pi/d pi/d],[point point],'--')
plot([-1*pi/d -0.9*pi/d],[point point_up])
plot([-1*pi/d -0.9*pi/d],[point point_down])
plot([0.9*pi/d pi/d],[point_up point])
plot([0.9*pi/d pi/d],[point_down point])
text(-0.4*pi/d,1.2,'Visible region')

plot([pi/d 3*pi/d],[point point],'--')
plot([pi/d 1.1*pi/d],[point point_up])
plot([pi/d 1.1*pi/d],[point point_down])
plot([2.9*pi/d 3*pi/d],[point_up point])
plot([2.9*pi/d 3*pi/d],[point_down point])
text(1.6*pi/d,1.2,'Virtual region')

subplot(4,1,2)
plot(psi,B);
hold on
axis([u1*pi u2*pi y1 y2])
set(gca,'XTick',[-pi 0 pi 2*pi 3*pi])
set(gca,'YTick',[0 0.5 1])
set(gca,'XTickLabelMode','manual')
%   set(gca,'XTickLabel',{'-pi/d';'0';'pi/d';'2pi/d';'3pi/d'})
c=0.18;
b=0.15;
text(-pi-c,y1-b,'-\pi')
text(0-0.1,y1-b,'0')
text(pi-0.1,y1-b,'\pi')
text(2*pi-c,y1-b,'2\pi')
text(3*pi-c,y1-b,'3\pi')

grid on
legend('{\it \psi}-space')

subplot(4,1,3)
plot(psi/pi,B);
hold on
axis([u1 u2 y1 y2])
set(gca,'XTick',[-1 0 1 2 3])
set(gca,'YTick',[0 0.5 1])
grid on
legend('{\it u}-space')
theta = [-180:0.1:360];
vv = exp(j*2*d*n*pi*cos(theta/180*pi));
B = (abs(w'*vv));   
subplot(4,1,4)
plot(theta,B);
hold on
axis([-180 180 y1 y2])
set(gca,'XTick',[-180 -90 0 90 180])
set(gca,'YTick',[0 0.5 1])
set(gca,'XTickLabelMode','manual')

%   set(gca,'XTickLabel',{'180';'90';'0';'-90';'-180'})
c=8;
b=0.15;
text(-180-c,y1-b,'180^o')
text(-90-c,y1-b,'90^o')
text(0-4,y1-b,'0^o')
text(90-c,y1-b,'-90^o')
text(180-c,y1-b,'-180^o')

grid on
legend('\theta-space')