%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 4.1.4
% K. Bell 10/19/99
% updated by K. Bell 11/3/00
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all


Nx = 11;
Ny = 11;
dx = 0.5;
dy = 0.5;
Dx = [-(Nx-1)/2:1:(Nx-1)/2];
Dy = [-(Ny-1)/2:1:(Ny-1)/2];
ux = [-1:0.025:1];
uy = [-1:0.025:1];

D = zeros(Nx*Ny,2);
for m=1:Ny
   D((m-1)*Nx+1:m*Nx,:) = [dx*Dx.' dy*Dy(m)*ones(Nx,1)];
end
C = zeros(Nx,Ny);
for k=1:Nx*Ny;
   for l=1:Nx*Ny
      dp = D(k,:)-D(l,:);
      C(k,l) = sinc(2*sqrt(dp*dp'));
   end
end

nx = length(ux);
ny = length(uy);

BP = zeros(nx,ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dolph-Chebychev
R = 10^(30/20);
x0 = cosh(acosh(R)/(Nx-1));
wt = poly(exp(j*2*acos(cos((2*[1:1:Nx-1]-1)*pi/(2*(Nx-1)))/x0))).';
wx = wt/sum(wt);

x0 = cosh(acosh(R)/(Ny-1));
wt = poly(exp(j*2*acos(cos((2*[1:1:Ny-1]-1)*pi/(2*(Ny-1)))/x0))).';
wy = wt.'/sum(wt);

w = reshape(wx*wy,Nx*Ny,1);
w = w/sum(w);

for n=1:nx
   for m=1:ny
      V = exp(j*2*pi*D*[ux(n);uy(m)]);
      BP(n,m) = w'*V;
   end
end

figure(1)
subplot(2,1,1)
G = 10*log10(abs(BP).^2);
G = max(G,-100);
mesh(ux,uy,G)
xlabel('ux')
ylabel('uy')
title('Prob. 4.1.4, separable Dolph-Chebychev weights')
axis([-1 1 -1 1 -100 0])

subplot(2,1,2)
contourf(ux,uy,G,[-10 -20 -30 -40 -50 -60 -70 -80 -90 -100])
axis('square')
colormap('gray')
colorbar
xlabel('ux')
ylabel('uy')
drawnow
set(gcf,'Paperposition',[0.25 1 8 9])

phi = 0;
% sin theta
u=[-1:0.001:1];
ux = u*cos(phi);
uy = u*sin(phi);
V = exp(j*2*pi*D*[ux;uy]);
Bc = w'*V;
figure(2)
subplot(3,1,1)
plot(u,10*log10(abs(Bc).^2))
axis([-1 1 -100 0])
grid on
title('Prob. 4.1.4, separable Dolph-Chebychev weights, \phi=0')
xlabel('ur')
ylabel('Beampattern (dB)')

phi = pi/4;
ux = u*cos(phi);
uy = u*sin(phi);
V = exp(j*2*pi*D*[ux;uy]);
Bc = w'*V;
subplot(3,1,2)
plot(u,10*log10(abs(Bc).^2))
axis([-1 1 -100 0])
grid on
title('\phi=\pi/4')
xlabel('ur')
ylabel('Beampattern (dB)')

phi = pi/2;
ux = u*cos(phi);
uy = u*sin(phi);
V = exp(j*2*pi*D*[ux;uy]);
Bc = w'*V;
subplot(3,1,3)
plot(u,10*log10(abs(Bc).^2))
axis([-1 1 -100 0])
title('\phi=\pi/2')
grid on
xlabel('ur')
ylabel('Beampattern (dB)')

set(gcf,'Paperposition',[0.25 1 8 9])

D = abs(w'*ones(Nx*Ny,1))^2/real(w'*C*w)
DI = 10*log10(D)

