%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.2.6
% K. Bell 11/17/99
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
ua = [0:-0.001:-0.1];
na = length(ua);
va = exp(j*pi*D*ua);

sigma_a = 10.^([-20 0 15 30]/10);
nSNR = length(sigma_a);
st = [':k';'--';'-.';'-k'];

%************************
% Source
%************************
us = 0.0;
AS = exp(j*pi*D*us);
sigma_n = 1;
SINRin = sigma_a/sigma_n;

wc = AS/N;

Tv = zeros(nSNR,na);
Av = zeros(nSNR,na);
Tp = zeros(nSNR,na);
Ap = zeros(nSNR,na);
for k=1:nSNR
   
   for n=1:na
      Rs = sigma_a(k)*va(:,n)*va(:,n)';
      Sn = sigma_n*eye(N);
      Sx = Rs+sigma_n*eye(N);
      Sninv = inv(Sn);
      Sxinv = inv(Sx);
      wv = Sninv*AS/real(AS'*Sninv*AS);
      wp = Sxinv*AS/real(AS'*Sxinv*AS);
      Tv(k,n) = real(wv'*wv);
      Av(k,n) = real((wv'*Rs*wv)/(wv'*Sn*wv))/SINRin(k);
      Tp(k,n) = real(wp'*wp);
      Ap(k,n) = real((wp'*Rs*wp)/(wp'*Sn*wp))/SINRin(k);
   end
   
end
figure(1)
subplot(2,1,1)
h1=plot(ua,10*log10(Ap(1,:)),'--');
hold on
h2=plot(ua,10*log10(Ap(2,:)),':');
h3=plot(ua,10*log10(Ap(3,:)),'-.');
h4=plot(ua,10*log10(Ap(4,:)),'-');
legend([h1 h2 h3 h4],'INR = -10 dB','INR = 0 dB','INR = 10 dB','INR = 20 dB')
xlabel('u_a')
ylabel('Array Gain (dB)')
title('Prob. 6.2.6, MPDR')
axis([-0.1 0 -40 10])
hold off
subplot(2,1,2)
h1=plot(ua,10*log10(Tp(1,:)),'--');
hold on
h2=plot(ua,10*log10(Tp(2,:)),':');
h3=plot(ua,10*log10(Tp(3,:)),'-.');
h4=plot(ua,10*log10(Tp(4,:)),'-');
legend([h1 h2 h3 h4],'INR = -10 dB','INR = 0 dB','INR = 10 dB','INR = 20 dB')
xlabel('u_a')
ylabel('Sensitivity (T) (dB)')
hold off
axis([-0.1 0 -15 30])
set(gcf,'Paperposition',[0.25 1 8 9])
figure(2)
subplot(2,1,1)
h1=plot(ua,10*log10(Av(1,:)),'--');
hold on
h2=plot(ua,10*log10(Av(2,:)),':');
h3=plot(ua,10*log10(Av(3,:)),'-.');
h4=plot(ua,10*log10(Av(4,:)),'-');
legend([h1 h2 h3 h4],'INR = -10 dB','INR = 0 dB','INR = 10 dB','INR = 20 dB')
xlabel('u_a')
ylabel('Array Gain (dB)')
title('Prob. 6.2.6, MVDR')
axis([-0.1 0 -40 10])
hold off
subplot(2,1,2)
h1=plot(ua,10*log10(Tv(1,:)),'--');
hold on
h2=plot(ua,10*log10(Tv(2,:)),':');
h3=plot(ua,10*log10(Tv(3,:)),'-.');
h4=plot(ua,10*log10(Tv(4,:)),'-');
legend([h1 h2 h3 h4],'INR = -10 dB','INR = 0 dB','INR = 10 dB','INR = 20 dB')
xlabel('u_a')
ylabel('Sensitivity (T) (dB)')
hold off
axis([-0.1 0 -15 30])
set(gcf,'Paperposition',[0.25 1 8 9])

