%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 4.1.4x
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

psix = [-(Nx-1)/2:1:(Nx-1)/2]*2*pi/Nx;
psiy = [-(Ny-1)/2:1:(Ny-1)/2]*2*pi/Ny;
nx = length(psix);
ny = length(psiy);

B = zeros(nx*ny,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=Nx-1;
R = 10^(30/20);
x0 = cosh(acosh(R)/m);

% Tseng-Cheng/Baklanov
V = zeros(Nx*Ny,Nx*Ny); 
B = zeros(Nx*Ny,1);
for n=1:nx
   for k=1:ny
      x = x0*cos(psix(n)/2)*cos(psiy(k)/2);  
      s= j*sqrt(1-x.^2);
      ind = (n-1)*ny+k;
      V(:,ind) = exp(j*D*[psix(n)/dx;psiy(k)/dx]);
      B(ind) = 0.5*((x+s).^m + (x-s).^m)/R;
   end
end
w1 = V'*B/(Nx*Ny);

% McClellan
V = zeros(Nx*Ny,Nx*Ny); 
B = zeros(Nx*Ny,1);
for n=1:nx
   for k=1:ny
      x = x0*cos(0.5*acos(cos(psix(n))*cos(psiy(k))));  
      s= j*sqrt(1-x.^2);
      ind = (n-1)*ny+k;
      V(:,ind) = exp(j*D*[psix(n)/dx;psiy(k)/dx]);
      B(ind) = 0.5*((x+s).^m + (x-s).^m)/R;
   end
end
w2 = V'*B/(Nx*Ny);


nx = length(ux);
ny = length(uy);

for q=1:2
   if q==1
      w = w1;
   else
      w = w2;
   end
   
   B2 = zeros(nx,ny);
   for n=1:nx
      for k=1:ny
         V = exp(j*2*pi*D*[ux(n);uy(k)]);
         B2(k,n)= w'*V;
      end
   end
   
   figure
   subplot(2,1,1)
   G = 10*log10(abs(B2).^2);
   G = max(G,-100);
   mesh(ux,uy,G)
   xlabel('ux')
   ylabel('uy')
   axis([-1 1 -1 1 -100 0])
   if q==1
      title('Prob. 4.1.4, Tseng-Cheng/Baklanov Dolph-Chebychev weights')
   else
      title('Prob. 4.1.4, McClellan Dolph-Chebychev weights')
   end
   
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
   uxx = u*cos(phi);
   uyy = u*sin(phi);
   V = exp(j*2*pi*D*[uxx;uyy]);
   Bc = w'*V;
   
   figure
   subplot(3,1,1)
   plot(u,10*log10(abs(Bc).^2))
   axis([-1 1 -100 0])
   grid on
   if q==1
      title('Prob. 4.1.4, Tseng-Cheng/Baklanov Dolph-Chebychev weights, \phi=0')
   else
      title('Prob. 4.1.4, McClellan Dolph-Chebychev weights, \phi=0')
   end
   xlabel('ur')
   ylabel('Beampattern (dB)')
   
   phi = pi/4;
   uxx = u*cos(phi);
   uyy = u*sin(phi);
   V = exp(j*2*pi*D*[uxx;uyy]);
   Bc = w'*V;
   subplot(3,1,2)
   plot(u,10*log10(abs(Bc).^2))
   axis([-1 1 -100 0])
   grid on
   title('\phi=\pi/4')
   xlabel('ur')
   ylabel('Beampattern (dB)')
   
   phi = pi/2;
   uxx = u*cos(phi);
   uyy = u*sin(phi);
   V = exp(j*2*pi*D*[uxx;uyy]);
   Bc = w'*V;
   subplot(3,1,3)
   plot(u,10*log10(abs(Bc).^2))
   axis([-1 1 -100 0])
   title('\phi=\pi/2')
   grid on
   xlabel('ur')
   ylabel('Beampattern (dB)')
   
   set(gcf,'Paperposition',[0.25 1 8 9])
   DIR = abs(w'*ones(Nx*Ny,1))^2/real(w'*C*w)
   DI = 10*log10(DIR)
   
end % for q
