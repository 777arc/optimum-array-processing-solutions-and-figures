%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 5.4.3
% K. Bell 11/17/00
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 9;
d = 0.5;

DC = [0 0;0 1;0 2;0 3;0 4;1 0;2 0;3 0;4 0]*d;
DS = [-1 1;0 1;1 1;-1 0;0 0;1 0;-1 -1; 0 -1;1 -1]*d;
w=ones(N,1)/N;

figure(1)
clf
plot(DC(:,1),DC(:,2),'*')
axis('square')
hold on
plot(DS(:,1),DS(:,2),'or')
hold off

CC = zeros(N,N);
for k=1:N;
   for l=1:N
      dp = DC(k,:)-DC(l,:);
      CC(k,l) = sinc(2*sqrt(dp*dp'));
   end
end
CS = zeros(N,N);
for k=1:N;
   for l=1:N
      dp = DS(k,:)-DS(l,:);
      CS(k,l) = sinc(2*sqrt(dp*dp'));
   end
end
ux = [-1:0.025:1];
uy = [-1:0.025:1];
nx = length(ux);
ny = length(uy);

BC = zeros(ny,nx);
BS = zeros(ny,nx);
for n=1:nx
   for k=1:ny
      VC = exp(j*2*pi*DC*[ux(n);uy(k)]);
      BC(k,n)= w'*VC;
      VS = exp(j*2*pi*DS*[ux(n);uy(k)]);
      BS(k,n)= w'*VS;
   end
end

for k=1:2
   if k==1
      D = DC;
      G = 10*log10(abs(BC).^2);
      figure(2)
   else
      D=DS;
      G = 10*log10(abs(BS).^2);
      
      figure(4)
   end
   ux = [-1:0.025:1];
   uy = [-1:0.025:1];
   
   subplot(2,1,1)
   G = max(G,-60);
   mesh(ux,uy,G)
   xlabel('ux')
   ylabel('uy')
   if k==1
      title('Prob. 5.4.3, Crossed Array, Uniform weighting')
   else
      title('Prob. 5.4.3, Square Array, Uniform weighting')
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
   
   
   u=[-1:0.001:1];
   phi = 0;
   ux = u*cos(phi);
   uy = u*sin(phi);
   V = exp(j*2*pi*D*[ux;uy]);
   Bc = w'*V;
   if k==1
      figure(3)
   else
      figure(5)
   end
   
   subplot(3,1,1)
   plot(u,10*log10(abs(Bc).^2))
   axis([-1 1 -60 0])
   grid on
   if k ==1
      title('Prob. 5.4.3, Crossed Array, Uniform weighting, \phi=0')
   else
      title('Prob. 5.4.3, Square Array, Uniform weighting, \phi=0')
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
   
   phi = -pi/4;
   ux = u*cos(phi);
   uy = u*sin(phi);
   V = exp(j*2*pi*D*[ux;uy]);
   Bc = w'*V;
   subplot(3,1,3)
   plot(u,10*log10(abs(Bc).^2))
   axis([-1 1 -60 0])
   title('\phi=-\pi/4')
   grid on
   xlabel('ur')
   ylabel('Beampattern (dB)')
   
   set(gcf,'Paperposition',[0.25 1 8 9])
   if k==1
      C=CC;
   else
      C=CS;
   end
  
   DIR = abs(w'*ones(N,1))^2/real(w'*C*w)
   DI = 10*log10(DIR)
end


