%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 5.5.25
% K. Bell 11/9/99
% updated by K. Bell 11/17/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%************************
% Array
%************************
N = 10;                                 % Elements in array
d= 0.5;                                % sensor spacing half wavelength wrt wc
D = d*[-(N-1)/2:1:(N-1)/2].';
u = [-1:0.005:1];
nx = length(u);
vv = exp(j*2*pi*D*u);

%************************
% Source
%************************

S2 = [0 10 20 30];

for n2 = 1:4
   sigma_s = [1 0;0 10.^(S2(n2)/10)];
   du = [[0.01:0.005:0.19] [0.2:0.025:0.9]];
   
   ns = length(du);
   us1 = -du/2;
   us2 = du/2;
   V1 = exp(j*2*pi*D*us1);
   V2 = exp(j*2*pi*D*us2);
   
   val1 = zeros(1,ns);
   val2 = zeros(1,ns);
   
   
   for n=1:ns
      V = [V1(:,n) V2(:,n)];
      R = V*sigma_s*V';
      [v,lam] = eig(R);
      [lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
      v       = v(:,ind);                      % arrange eigenvectors in same order
      v1 = v(:,N);
      v2 = v(:,N-1);
      val1(n) = real(lam(N))/N;
      val2(n) = real(lam(N-1))/N;
   end
   
   figure(1)
   semilogx(du,10*log10(val1),'-')
   hold on
   semilogx(du,10*log10(val2),'--g')
   q=axis;
   axis([q(1:2) -30 35])
   grid off
   xlabel('|\Delta u|')
   ylabel('Eigenvalues (dB)')
   title(['Prob. 5.5.25'])
   hold on
   
   du = [0.05 0.3 0.7];
   
   ns = length(du);
   us1 = -du/2;
   us2 = du/2;
   V1 = exp(j*2*pi*D*us1);
   V2 = exp(j*2*pi*D*us2);
   figure
   
   for n=1:ns
      V = [V1(:,n) V2(:,n)];
      R = V*sigma_s*V';
      [v,lam] = eig(R);
      [lam,ind] = sort(diag(abs(lam)));        % sort eigenvalues in ascending order 
      v       = v(:,ind);                      % arrange eigenvectors in same order
      v1 = v(:,N);
      v2 = v(:,N-1);
      B1 = v1'*vv;
      B2 = v2'*vv;
      subplot(3,1,n)
      h1=plot(u,10*log10(abs(B1).^2),'-');
      hold on
      h2=plot(u,10*log10(abs(B2).^2),'--');
      plot(us1(n)*[1 1],[-50 10],'-')
      plot(us2(n)*[1 1],[-50 10],'-')
      hold off
      axis([-1 1 -30 10])
      ylabel('Eigenbeams (dB)')
      if n==1
         title(['Prob. 5.5.25, \Delta u = ' num2str(du(n)) ', S2 = ' num2str(10*log10(sigma_s(2,2))) 'dB'])
      else
         title(['\Delta u = ' num2str(du(n))])
      end
   end
   xlabel('u')
   set(gcf,'Paperposition',[0.25 1 8 9])
end



