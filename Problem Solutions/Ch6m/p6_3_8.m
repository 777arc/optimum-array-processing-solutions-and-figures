%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.3.8
% K. Bell 11/19/00
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
ua = [[0.0:0.001:0.04] [0.05:0.01:0.25]]*BWNN;
na = length(ua);
nu=length(u);
vv = exp(j*pi*D*u);
va = exp(j*pi*D*ua);

sigma_s = 10.^([-10 0 10 30]/10);

%************************
% Source
%************************
us = 0.0;
AS = exp(j*pi*D*us);

%************************
% Interference
%************************
uI = [0.3 0.5];
vI = exp(j*pi*D*uI);
sigma_i = 10.^([20 20]/10);
sigma_n = 1;
Sn = vI*[sigma_i(1) 0;0 sigma_i(2)]*vI'+sigma_n*eye(N);
IN = trace(Sn)/N;

A = zeros(2,na);

for n=1:na
   for k=1:4
      Ss = sigma_s(k)*va(:,n)*va(:,n)';
      SINRin = sigma_s(k)/IN;

      Sx = Ss+Sn;
      Sxinv = inv(Sx);
      w = Sxinv*AS/real(AS'*Sxinv*AS);
      A(k,n) = 10*log10(real((w'*Ss*w)/(w'*Sn*w))/SINRin);
      
      if ua(n)==0.02*BWNN | ua(n)==0.05*BWNN |ua(n)==0.1*BWNN
         if k==1
            figure
            set(gcf,'Paperposition',[0.25 1 8 9])
         end
         B = 10*log10(abs(w'*vv).^2);

         subplot(4,1,k)
         plot(u,B,'-g');
         hold on
         plot(ua(n)*[1 1],[-60 30],'c')
         plot(uI(1)*[1 1],[-60 30],'r')
         plot(uI(2)*[1 1],[-60 30],'r')
         xlabel('u')
         ylabel('Beampattern (dB)')
         grid on
         hold off
         if k==1
            title(['Prob. 6.3.8, SNR = ' num2str(10*log10(sigma_s(k))) 'dB, ua = ' num2str(ua(n)/BWNN) ' BWNN'])
         else
            title(['SNR = ' num2str(10*log10(sigma_s(k))) 'dB, ua = ' num2str(ua(n)/BWNN) ' BWNN'])
         end
         axis([-1 1 -60 30])
         set(gca,'Ytick',[-60 -40 -20 0 20])
      end
   end
end

figure
h1=plot(ua/BWNN,A(1,:),'-');
hold on
h2=plot(ua/BWNN,A(2,:),'--');
h3=plot(ua/BWNN,A(3,:),':');
h4=plot(ua/BWNN,A(4,:),'-.');
hold off
title('Prob. 6.3.8')
legend('SNR = -10 dB','SNR = 0 dB','SNR = 10 dB','SNR = 20 dB')
axis([0 0.25 -20 50])
grid on
xlabel('u_a/BWNN')
ylabel('Array gain (dB)')

