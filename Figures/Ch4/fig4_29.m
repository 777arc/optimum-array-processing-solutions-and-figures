%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.29
% Beampattern with Riblet weighting, d=lambda/4, SLL=-20dB
% Xiaomin Lu 11-2-98	
% Updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
% Function called: cheby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%

clear all
close all

N = 11;
[ux,uy] = meshgrid(-1:1/40:1);
m = 1;
dB = 20;
d = 1/4;

R = 10^(dB/20);
x0 = cosh( 2/(N-1)*acosh(R));
alpha = 2*pi*d;

x = 1/(1-cos(alpha))*( (x0+1)*cos(alpha*ux).*cos(alpha*uy)-1-x0*cos(alpha));
beam = cheby((N-1)/2,x);
beam = abs(beam)/max(max(abs(beam)));
mesh(ux,uy,max(20*log10(beam),-80))
axis([-1 1 -1 1 -80 0])
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14);
%title('Riblet pattern using cos(\phi)=cos(\phi_x)*cos(\phi_y), N=11, d=\lambda/4')
figure

u = -1:1/100:1;
t = 1;
for phi = [0 30 60 90]
    ux = u*cos(phi/180*pi);
    uy = u*sin(phi/180*pi);
    x = 1/(1-cos(alpha))*( (x0+1)*cos(alpha*ux).*cos(alpha*uy)-1-x0*cos(alpha));
    beam = cheby((N-1)/2,x);
    beam = abs(beam)/max(abs(beam));
    table(t,:) = 20*log10(abs(beam));
    t = t+1;
end
subplot(2,2,1)
plot(u,table(1,:))
grid
axis([-1 1 -60 10])
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
text(-0.9,4,'\phi = 0^o','Fontsize',12)
%text(-1,18,'               Riblet Pattern cuts, N = 11')

subplot(2,2,2)
plot(u,table(2,:))
grid
axis([-1 1 -60 10])
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
text(-0.9,4,'\phi = 30^o','Fontsize',12)

subplot(2,2,3)
plot(u,table(3,:))
grid
axis([-1 1 -60 10])
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
text(-0.9,4,'\phi = 60^o','Fontsize',12)

subplot(2,2,4)
plot(u,table(4,:))
grid
axis([-1 1 -60 10])
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
text(-0.9,4,'\phi = 90^o','Fontsize',12)



