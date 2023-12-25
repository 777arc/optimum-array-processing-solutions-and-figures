%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.48
% Beampattern of sampled Taylor weighting
% Xiaomin Lu 11-2-98	
% Updated by K. Bell 10/11/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
% Function called: taylor_circ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%
clear all
close all

N = 20;
scale = 5;

i =1;
for t = [0 15/180*pi 30/180*pi 45/180*pi]
    [theta,phi] =meshgrid(0:pi/400:pi/2,t);
    
    ux = sin(theta).*cos(phi);
    uy = sin(theta).*sin(phi);
    
    
    Beam = 0;
    for n = 1:N/2
        for m = 1:N/2
            s1 = n-0.5;
            s2 = m-0.5;
            t1 = sqrt(s1^2+s2^2);
            if (t1<=N/2)
                Pmn = pi/2/scale*t1;
                w = taylor_circ(Pmn);
                Beam = Beam + 4*w*cos(s1*pi*ux).*cos(s2*pi*uy);
            end
        end
    end
    
    Beam = 20*log10( abs(Beam)/max(abs(Beam)) );
    table(i,:) = Beam;
    i = i+1;
end


xscale = 0:pi/400:pi/2;
subplot(2,2,1)
plot(sin(xscale),table(1,:))
axis([0 1 -50 0])
xlabel('\it u_r','Fontsize',14)
ylabel('Beam pattern(dB)','Fontsize',14)
grid
title('\phi=0^o','Fontsize',14);
%text(0,8,'Beam of sampled Taylor wght, N=20, SLL=-20dB, nBar=6, R=5*\lambda');


subplot(2,2,2)
plot(sin(xscale),table(2,:))
axis([0 1 -50 0])
xlabel('\it u_r','Fontsize',14)
ylabel('Beam pattern(dB)','Fontsize',14)
grid
title('\phi=15^o','Fontsize',14);

subplot(2,2,3)
plot(sin(xscale),table(3,:))
axis([0 1 -50 0])
xlabel('\it u_r','Fontsize',14)
ylabel('Beam pattern(dB)','Fontsize',14)
grid
title('\phi=30^o','Fontsize',14);

subplot(2,2,4)
plot(sin(xscale),table(4,:))
axis([0 1 -50 0])
xlabel('\it u_r','Fontsize',14)
ylabel('Beam pattern(dB)','Fontsize',14)
grid
title('\phi=45^o','Fontsize',14);
