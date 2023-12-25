%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 5.5.1
% K. Bell 11/9/99
% updated by K. Bell 11/20/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
%************************
% Array
%************************
Nx = 10;                                 % Elements in array
dx= 0.5;                                % sensor spacing half wavelength wrt wc
Ny = 10;                                 % Elements in array
dy= 0.5;                                % sensor spacing half wavelength wrt wc
Dx = [-(Nx-1)/2:1:(Nx-1)/2];
Dy = [-(Ny-1)/2:1:(Ny-1)/2];
ux = [-1:0.025:1];
uy = [-1:0.025:1];
nx = length(ux);
ny = length(uy);

N = Nx*Ny;

D = zeros(Nx*Ny,2);
for m=1:Ny
   D((m-1)*Nx+1:m*Nx,:) = [dx*Dx.' dy*Dy(m)*ones(Nx,1)];
end

sigma_s = 10.^([0 0]/10);
sigma_s = [sigma_s(1) 0;0 sigma_s(2)];

%************************
% Source
%************************
phi = [0 pi/4 pi/3 pi/2];
du = [[0.01:0.005:0.1] [0.2:0.025:1]];
dup = [0.05 0.3 0.7];

ns = length(du);
val1 = zeros(1,ns);
val2 = zeros(1,ns);
exact1 = zeros(1,ns);
exact2 = zeros(1,ns);
for p = 1:4
   ['phi = ' num2str(180*phi(p)/pi)]
   usx1 = -0.5*du*cos(phi(p));
   usy1 = -0.5*du*sin(phi(p));
   usx2 = 0.5*du*cos(phi(p));
   usy2 = 0.5*du*sin(phi(p));
   V1 = exp(j*2*pi*D*[usx1;usy1]);
   V2 = exp(j*2*pi*D*[usx2;usy2]);
   
   for n=1:ns
      R = [V1(:,n) V2(:,n)]*sigma_s*[V1(:,n) V2(:,n)]';
      [v,lam] = eig(R);
      [lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
      v       = v(:,ind);                      % arrange eigenvectors in same order
      v1 = v(:,N);
      v2 = v(:,N-1);
      val1(n) = real(lam(N))/N;
      val2(n) = real(lam(N-1))/N;
      B12 = V1(:,n)'*V2(:,n)/N;
      exact1(n) = 1 + abs(B12);
      exact2(n) = 1 - abs(B12);
      e1 = (V1(:,n)+sign(B12)*V2(:,n))/sqrt(2*N*(1+sign(B12)*real(B12)));
      e2 = (V1(:,n)-sign(B12)*V2(:,n))/sqrt(2*N*(1-sign(B12)*real(B12)));
   end
   
   figure(1)
   subplot(4,1,p)
   semilogx(du,10*log10(val1),'-')
   hold on
   semilogx(du,10*log10(val2),'--')
   semilogx(du,10*log10(exact1),'-r')
   semilogx(du,10*log10(exact2),'--r')
   hold off
   v=axis;
   axis([v(1:2) -30 10])
   grid on
   ylabel('Eigenvalues (dB)')
   if p==1
      title('Prob. 5.5.1, Sources on u_x-axis')
   elseif p==2
      title('Sources on u_y=u_x axis')
   elseif p==3
      title('Sources on u_y=sqrt(3)*u_x axis')
   else
      title('Sources on u_y-axis')
      xlabel('|\Delta u|')
   end
   
   BP1 = zeros(nx,ny);
   BP2 = zeros(nx,ny);
   usx1 = -0.5*dup*cos(phi(p));
   usy1 = -0.5*dup*sin(phi(p));
   usx2 = 0.5*dup*cos(phi(p));
   usy2 = 0.5*dup*sin(phi(p));
   V1 = exp(j*2*pi*D*[usx1;usy1]);
   V2 = exp(j*2*pi*D*[usx2;usy2]);
   
   for n=1:3
      R = [V1(:,n) V2(:,n)]*sigma_s*[V1(:,n) V2(:,n)]'+eye(N);
      [v,lam] = eig(R);
      [lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
      v       = v(:,ind);                      % arrange eigenvectors in same order
      v1 = v(:,N);
      v2 = v(:,N-1);
      for k=1:nx
         for m=1:ny
            VV = exp(j*2*pi*D*[ux(k);uy(m)]);
            BP1(m,k) = v1'*VV;
            BP2(m,k) = v2'*VV;
         end
      end
      G1 = 10*log10(abs(BP1).^2);
      G1 = max(G1,-100);
      G2 = 10*log10(abs(BP2).^2);
      G2 = max(G2,-60);
      
      figure(p+1)
      set(gcf,'Paperposition',[0.25 1 7.5 9])
      subplot(3,2,(n-1)*2 +1)
      contourf(ux,uy,G1,[0 -10 -20 -30 -40 -50 -60])
      axis('square')
      colormap('gray')
      colorbar
      xlabel('ux')
      ylabel('uy')
      title(['Eigenbeam 1, |\Delta u| = ' num2str(dup(n))])
      subplot(3,2,(n-1)*2+2)
      contourf(ux,uy,G2,[0 -10 -20 -30 -40 -50 -60])
      axis('square')
      colormap('gray')
      colorbar
      xlabel('ux')
      ylabel('uy')
      title(['Eigenbeam 2, |\Delta u| = ' num2str(dup(n))])
   end
   drawnow
end
figure(1)
set(gcf,'Paperposition',[0.25 1 8 9])


