%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.2.4
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
uI = [[0:0.001:0.1] [0.105:0.005:1]]*BWNN;
nI = length(uI);
vI = exp(j*pi*D*uI);

sigma_i = 10.^([-10 0 10 20]/10);
nINR = length(sigma_i);
st = [':k';'--';'-.';'-k'];

%************************
% Source
%************************
us = 0.0;
AS = exp(j*pi*D*us);
sigma_n = 1;
sigma_s = 1;
SINRin = sigma_s*ones(1,nINR)./(sigma_i+sigma_n*ones(1,nINR));

wc = AS/N;
Rs = sigma_s*AS*AS';

T = zeros(nINR,nI);
Tc1 = zeros(nINR,nI);
Tc2 = zeros(nINR,nI);
A = zeros(nINR,nI);
Ac = zeros(nINR,nI);
for k=1:nINR
   g1 = sigma_i(k)/sigma_n;
   
   for n=1:nI
      Sn = sigma_i(k)*vI(:,n)*vI(:,n)'+sigma_n*eye(N);
      rho1s = vI(:,n)'*AS/N;
      as = (1-abs(rho1s).^2);
      Sninv = inv(Sn);
      w = Sninv*AS/real(AS'*Sninv*AS);
      T(k,n) = real(w'*w);
      Tc1(k,n) = (((1+N*g1)*as/(1+N*g1*as))^2)/N;
      Tc2(k,n) = (1+2*N*g1*as+N*N*g1*g1*as)/(N*(1+N*g1*as)^2);
      A(k,n) = real((w'*Rs*w)/(w'*Sn*w))/SINRin(k);
      AC(k,n) = N*(1+g1)*(1+N*g1*as)/(1+N*g1);
   end
end
figure(1)
subplot(2,1,1)
h1=plot(uI/BWNN,10*log10(A(1,:)),'--');
hold on
h2=plot(uI/BWNN,10*log10(A(2,:)),':');
h3=plot(uI/BWNN,10*log10(A(3,:)),'-.');
h4=plot(uI/BWNN,10*log10(A(4,:)),'-');
drawnow
xlabel('uI/BWNN')
ylabel('Array Gain (dB)')
title('Prob. 6.2.4')
axis([0 1 0 40])
hold off
subplot(2,1,2)
h1=plot(uI/BWNN,10*log10(T(1,:)),'--');
hold on
h2=plot(uI/BWNN,10*log10(T(2,:)),':');
h3=plot(uI/BWNN,10*log10(T(3,:)),'-.');
h4=plot(uI/BWNN,10*log10(T(4,:)),'-');
xlabel('uI/BWNN')
ylabel('Sensitivity (T) (dB)')
hold off
axis([0 1 -20 20])
set(gcf,'Paperposition',[0.25 1 8 9])

