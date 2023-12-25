%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 5.5.23
% K. Bell 11/23/98
% updated by K. Bell 11/17/00
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

sigma_s = 10^(0/10);

%************************
% Source
%************************
us = 0.0;
AS = exp(j*pi*D*us);

s = [0:0.01:1];
ns = length(s);
val = zeros(N,ns);
for n=1:ns
   r = exp(-0.5*([0:1:N-1]*pi*s(n)).^2);
   c = exp(-0.5*([0:1:N-1]*pi*s(n)).^2);
   Bs = toeplitz(c,r);
   Q = Bs.*(AS*AS');
   R = sigma_s*Q;
   
   [v,lam] = eig(R);
   [lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
   v       = v(:,ind);                      % arrange eigenvectors in same order
   val(:,n) = lam/N;
end
figure(1)
clf
plot(s,10*log10(val));
xlabel('\sigma_o')
ylabel('Eigenvalues (dB)')
grid on
hold off
axis([0 1.0 -50 0])
title('Prob. 5.5.23')

%%%%%%%%%%%%%%%%%%%%%%%%%

s = [0 0.05 0.1 0.2];
ns = length(s);
val = zeros(N,ns);
st = ['-b';'-r';'-c';'-m';'-g';'-k';':b';':r';':c';':m'];
for n=1:ns
   r = exp(-0.5*([0:1:N-1]*pi*s(n)).^2);
   c = exp(-0.5*([0:1:N-1]*pi*s(n)).^2);
   Bs = toeplitz(c,r);
   Q = Bs.*(AS*AS');
   R = sigma_s*Q;
   [v,lam] = eig(R);
   [lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
   v       = v(:,ind);                      % arrange eigenvectors in same order
   val(:,n) = lam/N;
   B = v'*vv;
   
   figure(2)
   subplot(4,1,n)
   B2 = B(N,:);
   plot(u,10*log10(abs(B2)),st(1,:));
   hold on
   ylabel('Eigenbeams')
   grid on
   t = sum(lam);
   rt = lam(N);
   for k=2:10
      rt = rt+lam(N-k+1);
      if rt<0.99*t
         plot(u,10*log10(abs(B(N-k+1,:)).^2),st(k,:));
      end
   end
   
   hold off
   axis([-1 1 -50 20])
   title(['\sigma_o = ' num2str(s(n))])
   
end
xlabel('u')

set(gcf,'Paperposition',[0.25 1 8 9])

