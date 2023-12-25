%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.3.2
% K. Bell 11/21/99
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
uI = [[1:-0.01:0.05] 0.0433 [0.04:-0.01:0]];
nI = length(uI);
nu=length(u);
vv = exp(j*pi*D*u);
vI = exp(j*pi*D*uI);

sigma_i = 10.^([10 50]/10);

%************************
% Source
%************************
us = 0.0;
AS = exp(j*pi*D*us);

sigma_n = 1;
sigma_s = 1;
SINRin = sigma_s./(sigma_i+sigma_n);

SL = -30;
nbar = 3;
x0 = cosh(acosh(10^(-SL/20))/(N-1));
cc=cos((2*[1:1:floor((N-1)/2)]-1)*0.5*pi/(N-1)).';
udc = acos(cc/x0)/(pi*d);
uu = [1:1:(N/2)-1].'/(N*d);
u_nbar = udc(nbar);       % u_nbar Chebychev
sigma = nbar/(N*d*u_nbar);
udcmod = udc*sigma;           % modified Chebychev
uz = [udcmod(1:nbar-1);uu(nbar:(N/2)-1)];
ut = [uz;-uz;-1];
w = poly(exp(j*2*pi*d*ut)).';
wc=w/sum(w);

Rs = sigma_s*AS*AS';

T = zeros(2,nI);
Tc = zeros(2,nI);
Tcn = zeros(2,nI);
A = zeros(2,nI);
Ac = zeros(2,nI);
Acn = zeros(2,nI);
Bc = 10*log10(abs(wc'*vv).^2);

for n=1:nI
   for k=1:2
      Sn = sigma_i(k)*vI(:,n)*vI(:,n)'+sigma_n*eye(N);
      Sninv = inv(Sn);
      w = Sninv*AS/real(AS'*Sninv*AS);
      PIo = eye(N)-vI(:,n)*vI(:,n)'/N;
      wcn = PIo*wc;
      T(k,n) = real(w'*w);
      Tc(k,n) = real(wc'*wc);
      Tcn(k,n) = real(wcn'*wcn);
      A(k,n) = 10*log10(real((w'*Rs*w)/(w'*Sn*w))/SINRin(k));
      Ac(k,n) = 10*log10(real((wc'*Rs*wc)/(wc'*Sn*wc))/SINRin(k));
      Acn(k,n) = 10*log10(real((wcn'*Rs*wcn)/(wcn'*Sn*wcn))/SINRin(k));
      B(k,:) = 10*log10(abs(w'*vv).^2);
      Bcn(k,:) = 10*log10(abs(wcn'*vv).^2);
      
      
      if uI(n)==0.5 | uI(n)==0.3 | uI(n)==0.18 |uI(n)==0.09 | uI(n)==0.0433 |uI(n)==0.02
         if k==1
            figure
            set(gcf,'Paperposition',[0.25 1 8 9])
         end
         subplot(2,1,k)
         h1=plot(u,B(k,:),'-b');
         hold on
         h2=plot(u,Bcn(k,:),'--g');
         h3=plot(u,Bc,'-.r');
         plot(us*[1 1],[-60 30],'c')
         plot(uI(n)*[1 1],[-60 30],'r')
         xlabel('u')
         ylabel('Beampattern (dB)')
         grid on
         hold off
         axis([-1 1 -60 20])
         set(gca,'Ytick',[-60 -40 -20 0 20])
         
         if k==1
            title(['Prob. 6.3.2, INR = ' num2str(10*log10(sigma_i(k))) 'dB, uI = ' num2str(uI(n))])
         else
            title(['INR = ' num2str(10*log10(sigma_i(k))) 'dB, uI = ' num2str(uI(n))])
         end
         legend([h1 h2 h3],['MVDR, AG =' sprintf('%2.1f',A(k,n)) ' dB'],...
            ['Vill. w/null, AG =' sprintf('%2.1f',Acn(k,n)) ' dB'],...
            ['Vill., AG =' sprintf('%2.1f',Ac(k,n)) ' dB'],2)
         drawnow
      end
   end
end
figure
for k=1:2
   subplot(2,1,k)
   h1=plot(uI,A(k,:),'-b');
   hold on
   h2=plot(uI,Acn(k,:),'--g');
   h3=plot(uI,Ac(k,:),'-.r');
   if k==1
      title(['Prob. 6.3.2, INR = ' num2str(10*log10(sigma_i(k))) 'dB'])
      axis([0 1 0 30])
   else
      title(['INR = ' num2str(10*log10(sigma_i(k))) 'dB'])
      axis([0 1 0 70])
   end
   hold off
   grid on
   xlabel('uI')
   ylabel('Array gain (dB)')
   legend([h1 h2 h3],'MVDR','Vill. w/null','Vill.',4)
   set(gcf,'Paperposition',[0.25 1 8 9])
end


   
   