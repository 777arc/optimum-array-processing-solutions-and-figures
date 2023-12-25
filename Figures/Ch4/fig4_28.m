%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4.28
% Sampling on 2-dimension T-C, using ifft to get weight
% and get the reconstructed beampattern using fft 
% Xiaomin Lu 11-2-98	
% updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
% Function called: cheby
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*****************************%

clear all
close all

N = 10;
M = 10;
dx = 0.5;
dy = 0.5;

dB = 20;
R = 10^(dB/20);
x0 = cosh(1/(N-1)*acosh(R));

[ux,uy] = meshgrid(-1:1/40:1);


for k1 = 0:N-1
    for k2 = 0:N-1
        Bpsi(k1+1,k2+1) = cheby(N-1,x0*cos((k1-(N-1)/2)*pi/N)*cos((k2-(M-1)/2)*pi/M));
    end
end

for k1 = 0:N-1
    for k2 = 0:N-1
        t1 = k1-(N-1)/2;
        t2 = k2-(M-1)/2;
        B(k1+1,k2+1) = exp(-j*(N-1)/2*(t1*2*pi/N)-j*(M-1)/2*(t2*2*pi/M))*Bpsi(k1+1,k2+1)';
    end
end

w = ifft2(B);
w = real(w)-j*imag(w);

for n = 0:N-1
    for m = 0:N-1
        w(n+1,m+1) = w(n+1,m+1)*exp(j*pi/N*(N-1)*(n+m));
    end
end 

w = abs(w)/max(max(abs(w)));
w(N/2+1:N,M/2+1:M)


Beam = 0;
for n = 0:N-1
    for m = 0:N-1
        Beam = Beam + w(n+1,m+1)*exp(j*n*2*pi*dx*ux+j*m*2*pi*dy*uy);
    end
end
Beam = Beam.*exp(-j*(N-1)/2*pi*ux-j*(N-1)/2*pi*uy);

Beam = 20*log10(abs(Beam)/max(max(abs(Beam))));
for k1 = 1:size(Beam,1)
    for k2 = 1:size(Beam,2)
        if (Beam(k1,k2)<-60) Beam(k1,k2)=-60; end
    end
end


mesh(ux,uy,Beam)
axis([-1 1 -1 1 -60 0])
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14);
%title('Regenerated D-C pattern, using IFFT, N=M=10, SLL=-20dB dx=dy=\lambda/2') 
figure

u = -1:1/100:1;
t = 1;
for phi = [0 30 60 90]
    ux = u*cos(phi/180*pi);
    uy = u*sin(phi/180*pi);
    
    Beam = 0;
    for n = 0:N-1
        for m = 0:N-1
            Beam = Beam + w(n+1,m+1)*exp(j*n*pi*ux+j*m*pi*uy);
        end
    end
    Beam = Beam.*exp(-j*(N-1)/2*pi*ux-j*(N-1)/2*pi*uy);
    table(t,:) = 20*log10(abs(Beam)/max(abs(Beam)));
    
    t = t+1;
end
subplot(2,2,1)
plot(u,table(1,:))
grid
axis([-1 1 -60 10])
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);
text(-0.9,4,'\phi = 0^o','Fontsize',12)
%text(-1,18,'                          Pattern cuts at several phis')

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




         
