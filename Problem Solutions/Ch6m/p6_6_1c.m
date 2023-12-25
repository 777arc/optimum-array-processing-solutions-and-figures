%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 6.6.1(c)
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

u1 = [0.0433 0.25 0.5]*BWNN;
n1=length(u1);
Sf = 10.^([-10:1:30]/10);
nf = length(Sf);
EA=zeros(n1,nf);

for q=1:n1
   us = [-1:0.01:1]*u1(q);
   du = 0.01*u1(q);
   AS = exp(j*2*pi*d*D*us);
   ns = length(us);
   for k=1:nf
      A = zeros(1,ns);
      for n=1:ns
         Sx = Sf(k)*AS(:,n)*AS(:,n)'+eye(N);
         Sxinv = inv(Sx);
         w = Sxinv*A0/(A0'*Sxinv*A0);
         A(n) = real(w'*AS(:,n)*AS(:,n)'*w)/real(w'*w);
      end
      EA(q,k) = sum(A)*du*0.5/u1(q);
   end
end

figure(1)
h1=plot(10*log10(Sf),10*log10(abs(EA(1,:))),'-');
hold on
h2=plot(10*log10(Sf),10*log10(abs(EA(2,:))),'--');
h3=plot(10*log10(Sf),10*log10(abs(EA(3,:))),'-.');
hold off
ylabel('E\{Array Gain\} (dB)')
xlabel('SNR (dB)')
axis([-10 30 -20 10])
grid on
legend([h1 h2 h3],['u1=' num2str(u1(1)/BWNN) ' BWNN'],['u1=' num2str(u1(2)/BWNN) ' BWNN'],['u1=' num2str(u1(3)/BWNN) ' BWNN']);
title('Problem 6.6.1(c)')