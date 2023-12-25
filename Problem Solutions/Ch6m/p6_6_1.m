%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.6.1
% K. Bell 11/30/98
% updated by K. Bell 11/26/00
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


%************************
% Source
%************************
u0 = 0;
A0 = exp(j*2*pi*d*D*u0);

us = [[0.0:0.01:0.49] [0.491:0.001:0.5]]*BWNN;
ns=length(us);
AS = exp(j*2*pi*d*D*us);
Sf = 10.^([-10:10:30]/10);
nf = length(Sf);
A=zeros(nf,ns);

for k=1:nf
   for n=1:ns
      Sx = Sf(k)*AS(:,n)*AS(:,n)'+eye(N);
      Sxinv = inv(Sx);
      w = Sxinv*A0/(A0'*Sxinv*A0);
      A(k,n) = real(w'*AS(:,n)*AS(:,n)'*w)/real(w'*w);
   end
end
figure(1)
set(gcf,'Paperposition',[0.25 1 8 9]);
subplot(2,1,1)
h1=plot(us/BWNN,10*log10(abs(A(1,:))),'-');
hold on
h2=plot(us/BWNN,10*log10(abs(A(2,:))),'--');
h3=plot(us/BWNN,10*log10(abs(A(3,:))),'-.');
h4=plot(us/BWNN,10*log10(abs(A(4,:))),':');
h5=plot(us/BWNN,10*log10(abs(A(5,:))),'-*');
hold off
ylabel('Array Gain (dB)')
title(['Problem 6.6.1'])
xlabel('ua/BWNN')
axis([0 0.5 -100 10])
grid on
legend([h1 h2 h3 h4 h5],['SNR=' num2str(10*log10(Sf(1))) ' dB'],['SNR=' num2str(10*log10(Sf(2))) ' dB'],['SNR=' num2str(10*log10(Sf(3))) ' dB'],['SNR=' num2str(10*log10(Sf(4))) ' dB'],['SNR=' num2str(10*log10(Sf(5))) ' dB'])

us = [0.25 0.0433]*BWNN;
ns=length(us);
AS = exp(j*2*pi*d*D*us);
Sf = 10.^([-10:1:30]/10);
nf = length(Sf);
A=zeros(ns,nf);

for k=1:nf
   for n=1:ns
      Sx = Sf(k)*AS(:,n)*AS(:,n)'+eye(N);
      Sxinv = inv(Sx);
      w = Sxinv*A0/(A0'*Sxinv*A0);
      A(n,k) = real(w'*AS(:,n)*AS(:,n)'*w)/real(w'*w);
   end
end
figure(1)
subplot(2,1,2)
h1=plot(10*log10(Sf),10*log10(abs(A(1,:))),'-');
hold on
h2=plot(10*log10(Sf),10*log10(abs(A(2,:))),'--');
hold off
ylabel('Array Gain (dB)')
xlabel('SNR (dB)')
axis([-10 30 -80 10])
grid on
legend([h1 h2],['ua=' num2str(us(1)/BWNN) ' BWNN'],['ua=' num2str(us(2)/BWNN) ' BWNN']);