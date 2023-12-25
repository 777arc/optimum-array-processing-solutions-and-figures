%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 4.4.1
% K. Bell 10/19/99
% updated by K. Bell 11/3/00
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

Nx = 9;NH = 61;
% kk = 0, 1, or 2
kk = 1;

dx = 0.5;
dy = 0.25*sqrt(3);
R = Nx*dx/2;

D = zeros(NH,2);
w=zeros(NH,1);
ind = 1;
for m=-(Nx-1)/2:(Nx-1)/2
   km = Nx-abs(m);
   for n= -(km-1)/2:(km-1)/2
      D(ind,:) = [dx*n dy*m];
      r = sqrt((dx*n)^2 +(dy*m)^2);
      w(ind) = (1-(r/R)^2)^kk;
      ind=ind+1;
   end
end

w=w/sum(w);

ux = [-1:0.025:1];
uy = [-1:0.025:1];
nx = length(ux);
ny = length(uy);

B = zeros(ny,nx);
for n=1:nx
   for k=1:ny
      V = exp(j*2*pi*D*[ux(n);uy(k)]);
      B(k,n)= w'*V;
   end
end

figure(1)
subplot(2,1,1)
G = 10*log10(abs(B).^2);
G = max(G,-60);
mesh(ux,uy,G)
xlabel('ux')
ylabel('uy')
if NH == 61 & kk==0
   title('Prob. 4.4.1(a), N=61, Uniform weighting')
elseif NH==61 & kk==1
   title('Prob. 4.4.1(b), N=61, Radial Taper')
else
   title('Prob. 4.4.1(c), N=61, Radial Taper Squared')
end

axis([-1 1 -1 1 -60 0])

subplot(2,1,2)
contourf(ux,uy,G,[-10 -20 -30 -40 -50 -60 -70 -80])
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
axis([-1 1 -60 0])
grid on
if NH == 61 & kk==0
   title('Prob. 4.4.1(a), N=61, Uniform weighting, \phi=0')
elseif NH==61 & kk==1
   title('Prob. 4.4.1(b), N=61, Radial Taper, \phi=0')
else
   title('Prob. 4.4.1(c), N=61, Radial Taper Squared, \phi=0')
end
xlabel('ur')
ylabel('Beampattern (dB)')

phi = pi/4;
ux = u*cos(phi);
uy = u*sin(phi);
V = exp(j*2*pi*D*[ux;uy]);
Bc = w'*V;
subplot(3,1,2)
plot(u,10*log10(abs(Bc).^2))
axis([-1 1 -60 0])
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
axis([-1 1 -60 0])
title('\phi=\pi/2')
grid on
xlabel('ur')
ylabel('Beampattern (dB)')

set(gcf,'Paperposition',[0.25 1 8 9])

