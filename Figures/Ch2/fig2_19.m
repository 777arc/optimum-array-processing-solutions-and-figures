%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.19
% Main lobe of beam pattern
% Kristine Bell
% Last updated 6/4/01, 10/4/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bp_sect - plots conventional bp with sector

clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uniform Linear Array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 10;                                % Elements in array
d = 0.5;                               % spacing wrt wavelength
beamwidth = 2/(N*d);                   % null-to-null (LL BW is half this)
D=d*[-(N-1)/2:1:(N-1)/2];              % element locations
u   = [-0.45:0.01:0.45]; 
AS  = exp(-j*2*pi*D'*0);               % BP points to 0
Au  = exp(-j*2*pi*D'*u);
B   = real(AS'*Au)/N;                          % BP

h=plot(u,B,'-');
set(h,'LineWidth',1.5)
hold on
% HPBW
I=find(B>=0.707);
plot([u(min(I)-1) u(max(I)+1)],[B(min(I)-1) B(max(I)+1)],'-')
plot(u(min(I)-1)*[1 1],B(min(I)-1)*[1 1]+0.03*[-1 1],'-')
plot(u(max(I)+1)*[1 1],B(max(I)+1)*[1 1]+0.03*[-1 1],'-')
plot([0.03 0.13],[0.72 0.845],'-')
h=text(0.155,0.84,'D');
set(h,'FontName','Symbol')
text(0.17,0.84,'{\itu}   = HPBW')
h=text(0.185,0.82,'1');
set(h,'Fontsize',10)
text(0.12, 0.7,'0.707')
% BW-NN
plot([-0.2 0.2],[-0.3 -0.3],'-')
plot(-0.2*[1 1],-0.3*[1 1]+0.03*[-1 1],'-')
plot(0.2*[1 1],-0.3*[1 1]+0.03*[-1 1],'-')
h=text(-0.05,-0.35,'D');
set(h,'FontName','Symbol')
text(-0.035,-0.35,'{\itu   = BW}')
h=text(-0.02,-0.37,'2');
set(h,'Fontsize',10)
h=text(0.04,-0.37,'\itNN');
set(h,'Fontsize',10)
% axes
plot([-1 1 ],[0 0],'-')
plot([0 0],[-0 1.1],'-')
% tick marks
plot(-0.4*[1 1],0.03*[-1 1])
plot(-0.2*[1 1],0.03*[-1 1])
plot(0.2*[1 1],0.03*[-1 1])
plot(0.4*[1 1],0.03*[-1 1])
text(-0.005,-0.05,'0')
%tick labels
xx=0.2;
yy=-0.1;
h=text(xx-0.01,yy+0.02,'l');
set(h,'FontName','Symbol')
plot([xx-0.02 xx+0.01],[yy-0.01 yy-0.01],'-')
h=text(xx-0.02,yy-0.05,'\itNd');

xx=-0.18;
yy=-0.1;
h=text(xx-0.01,yy+0.02,'l');
set(h,'FontName','Symbol')
plot([xx-0.02 xx+0.01],[yy-0.01 yy-0.01],'-')
h=text(xx-0.02,yy-0.05,'\itNd');
text(xx-0.04, yy-0.005,'-')

xx=0.41;
yy=-0.1;
h=text(xx-0.01,yy+0.02,'l');
set(h,'FontName','Symbol')
plot([xx-0.02 xx+0.01],[yy-0.01 yy-0.01],'-')
h=text(xx-0.02,yy-0.05,'\itNd');
text(xx-0.035, yy-0.005,'2')

xx=-0.4;
yy=-0.1;
h=text(xx-0.01,yy+0.02,'l');
set(h,'FontName','Symbol')
plot([xx-0.02 xx+0.01],[yy-0.01 yy-0.01],'-')
h=text(xx-0.02,yy-0.05,'\itNd');
text(xx-0.05, yy-0.005,'-2')

text(0.05,1.1,'\itB(u)')
hold off
axis([-0.42 0.42 -0.4 1.2])

%set(gca,'Box','off')
set(gca,'Visible','off')
