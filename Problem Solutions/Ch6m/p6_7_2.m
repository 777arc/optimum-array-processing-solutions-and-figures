%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.7.2
% K. Bell 11/30/98
% updated by K. Bell 11/29/00
% function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%************************
% Array
%************************
N = 10;                                 % Elements in array
d = 0.5;                                % sensor spacing half wavelength wrt wc
D = [-(N-1)/2:1:(N-1)/2].';
BWNN = 2/(N*d);
u=[-1:0.001:1];
nu=length(u);
vv = exp(j*pi*D*u);
uI = 0.3;
vI = exp(j*pi*D*uI);

sigma_i = 10.^([0 10 20]/10);
sigma_n = 1;
sigma_s = 10.^([-10 0 10 20]/10);
nf = length(sigma_s);
ni = length(sigma_i);

%DL = 10^(20/10);
DL=0;
%%%%%%%%%%%%%%%%%%%%
%(a) Directional
%uc = [-0.0433 0 0.0433];
%uc = [-0.0866 0 0.0866];
%uc = [-0.1 0 0.1];
%C = exp(j*pi*D*uc);
%Bcc = C'*ones(N,1)/N;
%g = [1 1 1].';
%g = Bcc;

%%%%%%%%%%%%%%%%%%%%%%%
% (b)Derivative
%uc = [0];
%C = [exp(j*pi*D*uc) j*pi*D.*exp(j*pi*D*uc) ((j*pi*D).^2).*exp(j*pi*D*uc)];
%Bcc = C'*ones(N,1)/N;
%g = [1;0;0];
%g = Bcc;

%%%%%%%%%%%%%%%%%%%%%%%
% (c) eigenvector
u_d = 0.1;
%u_d = 0.0866;
r = 2*u_d*sinc(-[0:1:N-1]*u_d);
c = 2*u_d*sinc([0:1:N-1]*u_d);
Q = toeplitz(c,r);
[v,lam] = eig(Q);
[lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
v       = v(:,ind);                      % arrange eigenvectors in same order

%p = 2*u_d*sinc(D*u_d);
p = Q*ones(N,1)/N;

C = [v(:,N) v(:,N-1) v(:,N-2)];          % first three eigenvectors
g = [p'*v(:,N)/lam(N) p'*v(:,N-1)/lam(N-1) p'*v(:,N-2)/lam(N-2)].';


%************************
% Source
%************************
us = [-0.25:0.005:0.25]*BWNN;
ns = length(us);
AS = exp(j*pi*D*us);

for q = 1:ni
   SINRo = zeros(nf,ns);
   for k=1:nf
      for n=1:ns
         Rs = sigma_s(k)*AS(:,n)*AS(:,n)';
         Snc = sigma_i(q)*vI*vI';
         Sn = Snc+sigma_n*eye(N);
         Sx = Rs+Sn;
         Sxinv = inv(Sx);
         w = Sxinv*C*inv(C'*Sxinv*C)*g;
         SINRo(k,n) = 10*log10(real(w'*Rs*w)/real(w'*Sn*w));
      end
   end
   
   figure(1)
   set(gcf,'Paperposition',[0.25 1 8 9]);
   subplot(ni,1,q)
   h1=plot(us/BWNN,SINRo(1,:),'-');
   hold on
   h2=plot(us/BWNN,SINRo(2,:),'--');
   h3=plot(us/BWNN,SINRo(3,:),'-.');
   h4=plot(us/BWNN,SINRo(4,:),':');
   if q==ni
      xlabel('u_a/BWNN')
         legend([h1 h2 h3 h4],['SNR = ' num2str(10*log10(sigma_s(1))) ' dB'],...
      ['SNR = ' num2str(10*log10(sigma_s(2))) ' dB'],...
      ['SNR = ' num2str(10*log10(sigma_s(3))) ' dB'],...
      ['SNR = ' num2str(10*log10(sigma_s(4))) ' dB']);

   end
   
   ylabel('SINRo (dB)')
   title(['INR = ' num2str(10*log10(sigma_i(q))) ' dB'])
   hold off
   axis([-0.25 0.25 -10 30])
end

%%%%%%%%%%
%Beampatterns
us = [-0.1 -0.05 0 0.05 0.1];
ns = length(us);
AS = exp(j*pi*D*us);
sigma_i = 10^(20/10);

%for q = 1:ni
figure(2)
set(gcf,'Paperposition',[0.25 1 8 9]);
for n=1:ns
   B = zeros(nf,nu);
   for k=1:nf
      Rs = sigma_s(k)*AS(:,n)*AS(:,n)';
      Snc = sigma_i*vI*vI';
      Sn = Snc+sigma_n*eye(N);
      Sx = Rs+Sn;
      Sxinv = inv(Sx);
      w = Sxinv*C*inv(C'*Sxinv*C)*g;
      
      B(k,:) = w'*vv;
   end
   subplot(3,2,n)
   h1=plot(u,10*log10(abs(B(1,:)).^2),'-');
   hold on
   h2=plot(u,10*log10(abs(B(2,:)).^2),'--');
   h3=plot(u,10*log10(abs(B(3,:)).^2),'-.');
   h4=plot(u,10*log10(abs(B(4,:)).^2),':');
   xlabel('u')
   
   if n==ns
      legend([h1 h2 h3 h4 ],['SNR = ' num2str(10*log10(sigma_s(1))) ' dB'],...
         ['SNR = ' num2str(10*log10(sigma_s(2))) ' dB'],...
         ['SNR = ' num2str(10*log10(sigma_s(3))) ' dB'],...
         ['SNR = ' num2str(10*log10(sigma_s(4))) ' dB']);
      
   end
   plot(uI*[1 1],[-60 20],'-')
   plot(us(n)*[1 1],[-60 20],'-')
   ylabel('Beampatterns (dB)')
   title(['INR = ' num2str(10*log10(sigma_i)) ' dB, u_a = ' num2str(us(n))])
   hold off
   axis([-1 1 -60 20])
end
   
%end

