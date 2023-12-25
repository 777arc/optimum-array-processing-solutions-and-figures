%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.24
% K. Bell 3/10/01, 11/12/01
% Function called: sinc
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
vv = exp(j*2*pi*d*D*u);

%************************
% Source
%************************

INR = 10^(20/10);
alpha = 1;

theta_s = [0:0.01:1]*pi;
us = cos(theta_s);
AS = exp(j*2*pi*d*D*us);
ns = length(us);
A = zeros(1,ns);
Ac = zeros(1,ns);
for n=1:ns
   p = [0:1:N-1];
   pI = [1:1:N-1];
   r = sinc(p*2*d)+[0 ((j*alpha)./(pi*pI*2*d)).*(sinc(pI*2*d)-cos(pi*2*d*pI))];
   Sn = INR*toeplitz(r,conj(r))+eye(N);
   Sninv = inv(Sn);
   Ac(n) =N*N*Sn(1,1)/ real(AS(:,n)'*Sn*AS(:,n)); 
   A(n) = real(AS(:,n)'*Sninv*AS(:,n))*Sn(1,1); 
end
figure
h1=plot(theta_s*180/pi,10*log10(A),'-');
hold on
h2=plot(theta_s*180/pi,10*log10(Ac),'--');
h=legend('MVDR', 'Conventional',4) ;
set(h,'Fontsize',12)
xlabel('\theta_{\its} (degrees)','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
grid on
hold off
%axis([0 180 5 20])
