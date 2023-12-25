%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 4.54, 4.55, 4.56
% Beampattern of 91-element uniform hexagonal array
% Xiaomin Lu 11-2-98
% Updated by K. Bell 10/11/00
% Lillian Xu updated 12/06/2000
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Beampattern for 91-element array, uniform weight

clear all
close all

Nx = 11;

[ux,uy] = meshgrid(-1:1/50:1);

Beam = 0;
for m = -(Nx-1)/2:(Nx-1)/2
    for n = -(Nx-1-abs(m))/2:(Nx-1-abs(m))/2
        Beam = Beam+exp(j*pi*n*ux+j*pi*m*sqrt(3)/2*uy);
    end
end
Beam = abs(Beam)/max(max(abs(Beam)));
Beam = 20*log10(Beam);
for k1 = 1:size(Beam,1)
    for k2 = 1:size(Beam,2)
        if(Beam(k1,k2)<-80) Beam(k1,k2)=-80; end
    end
end

mesh(ux,uy,Beam);
grid on
axis([-1 1 -1 1 -80 0])
%set(gca,'fontname','roman')
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14)
%title('Response of 91-element uniform Hexagonal array')
figure

cs = contour(ux,uy,Beam);
clabel(cs,'manual');
axis square
%set(gca,'fontname','roman')
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
%title('Contour plot of 91-element Hex Array, uniform weighting') 


u = -1:1/100:1;
k = 1;
for phi = [0 10 20 30]
    ux = u*cos(phi/180*pi);
    uy = u*sin(phi/180*pi);
    Beam = 0;
    for m = -(Nx-1)/2:(Nx-1)/2
        for n = -(Nx-1-abs(m))/2:(Nx-1-abs(m))/2
            Beam = Beam+exp(j*pi*n*ux+j*pi*m*sqrt(3)/2*uy);
        end
    end
    Beam = 20*log10(abs(Beam)/max(abs(Beam)));
    table(k,:) = Beam;
    k = k+1;
end
figure
subplot(2,2,1)
plot(u,table(1,:));
grid
axis([-1 1 -60 0])
%set(gca,'fontname','roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=0^o','Fontsize',14)
%text(-0.07,-72.5,'(a)')

subplot(2,2,2)
%figure
plot(u,table(2,:));
grid
axis([-1 1 -60 0])
%set(gca,'fontname','roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=10^o','Fontsize',14)
%text(-0.07,-72.5,'(b)')

subplot(2,2,3)
%figure
plot(u,table(3,:));
grid
axis([-1 1 -60 0])
%set(gca,'fontname','roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=20^o','Fontsize',14)
%text(-0.07,-72.5,'(c)')
subplot(2,2,4)
%figure
plot(u,table(4,:));
grid
axis([-1 1 -60 0])
%set(gca,'fontname','roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=30^o','Fontsize',14)
%text(-0.07,-72.5,'(d)')



