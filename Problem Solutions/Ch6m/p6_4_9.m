%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.4.9
% K. Bell 11/21/99
% updated by K. Bell 12/04/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%************************
% Array
%************************
N = 10;                                 % Elements in array
d = 0.5;                                % sensor spacing half wavelength wrt wc
D = [-(N-1)/2:1:(N-1)/2].';
BWNN = 2/(N*d);
u=[-1:0.001:1];
nu=length(u);
vv = exp(j*2*pi*d*D*u);


%************************
% Source
%************************
us = 0;
AS = exp(j*2*pi*d*D*us);
sigma_s = 1;
sigma_n = 1;

Rs = sigma_s*AS*AS';

uI = [0.3 0.2 0.5];
vI = exp(j*2*pi*d*D*uI);
uD = [0.1 0.1 0.2];
INR = 10.^([20 20 40]/10);
nI = length(INR);

   for q=1:nI
      r = sinc([0:1:N-1]*2*d*uD(q));
      Sn = INR(q)*toeplitz(r,r).*(vI(:,q)*vI(:,q)')+eye(N);
      Sninv = inv(Sn);
      w = Sninv*AS/real(AS'*Sninv*AS);
      %[trace(Sn)/N sigma_n+INR(q)]
      SINRin = sigma_s/(sigma_n+INR(q));
      AG = 10*log10(real((w'*Rs*w)/(w'*Sn*w))/SINRin);

      
      
      figure(1)
      set(gcf,'Paperposition',[0.25 1 8 9])
      subplot(nI,1,q)
      B = w'*vv;
      plot(u,10*log10(abs(B)));
      hold on
      plot([uI(q)-uD(q) uI(q)-uD(q) uI(q)+uD(q) uI(q)+uD(q)],[-40 0 0 -40],'r')
      plot([us us],[-40 10],'g')
      text(-0.95,5,['AG =' sprintf('%2.1f',AG) ' dB'])
      if q== nI
         xlabel('u')
      end
     
      ylabel('MVDR Beampattern (dB)')
      if q==1
         title(['Problem 6.4.9, u_I = ' num2str(uI(q)) ', u_{\Delta} = ' num2str(uD(q)) ', INR = ' num2str(10*log10(INR(q))) ' dB'])
      else
         title(['u_I = ' num2str(uI(q)) ', u_{\Delta} = ' num2str(uD(q)) ', INR = ' num2str(10*log10(INR(q))) ' dB'])
      end
      
      grid on
      axis([-1 1 -40 10])
      hold off
   end
set(gcf,'Paperposition',[0.25 1 8 9])
