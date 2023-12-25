%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures 4.57 & 4.58
% Beampattern for 91-element array with radial taper weighting, R = 2.75*lambda
% Xiaomin Lu  11-2-98
% Updated by K. Bell 10/11/00
% Updated by Lillian Xu 12/06/2000, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

Nx = 11;
R = 2.75;
[ux,uy] = meshgrid(-1:1/30:1);

Beam = 0;
for m = -(Nx-1)/2:(Nx-1)/2
    for n = -(Nx-1-abs(m))/2:(Nx-1-abs(m))/2
        w = 1-(4*n^2+3*m^2)/R^2/16;      
        Beam = Beam+w*exp(j*pi*n*ux+j*pi*m*sqrt(3)/2*uy);
    end
end
Beam = abs(Beam)/max(max(abs(Beam)));
Beam = 20*log10(Beam);

for k1 = 1:size(Beam,1)
    for k2 = 1:size(Beam,2)
        if(Beam(k1,k2)<-100) Beam(k1,k2)=-100; end
    end
end
figure
mesh(ux,uy,Beam);
grid on
axis([-1 1 -1 1 -100 0])
%et(gca, 'fontname', 'roman')
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14)
%title('Response of 91-ele Hex array, with 1st Radial Taper Wght, R=2.75\lambda')
figure


u = -1:1/100:1;
k = 1;
for phi = [0 10 20 30]
    ux = u*cos(phi/180*pi);
    uy = u*sin(phi/180*pi);
    Beam = 0;
    for m = -(Nx-1)/2:(Nx-1)/2
        for n = -(Nx-1-abs(m))/2:(Nx-1-abs(m))/2
            w = 1-(4*n^2+3*m^2)/R^2/16;      
            Beam = Beam+w*exp(j*pi*n*ux+j*pi*m*sqrt(3)/2*uy);
        end
    end
    Beam = 20*log10(abs(Beam)/max(abs(Beam)));
    table(k,:) = Beam;
    k = k+1;
end
subplot(2,2,1)
plot(u,table(1,:));
grid
axis([-1 1 -60 0])
%set(gca, 'fontname', 'roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=0^o','Fontsize',14)


subplot(2,2,2)
plot(u,table(2,:));
grid
axis([-1 1 -60 0])
%set(gca, 'fontname', 'roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=10^o','Fontsize',14)

subplot(2,2,3)
plot(u,table(3,:));
grid
axis([-1 1 -60 0])
%set(gca, 'fontname', 'roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=20^o','Fontsize',14)

subplot(2,2,4)
plot(u,table(4,:));
grid
axis([-1 1 -60 0])
%set(gca, 'fontname', 'roman')
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
title('\phi=30^o','Fontsize',14)




