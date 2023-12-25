%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.18
% Zhi Tian 7/12/99	
% Updated by K. Bell 10/13/00, 7/24/01, 10/23/01
% function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=10;
amf=[-(N-1)/2:(N-1)/2]';

u=[0.001:0.001:1];
SNRdb=[0:10:30];
ASNR=10.^(SNRdb/10)*N;

for i=1:length(ASNR)
   a=N*(1-2/ASNR(i));
   c1=[1 zeros(1,N-2) -a 0 a zeros(1,N-2) -1];
   c2=[1 zeros(1,N-2) a 0 -a zeros(1,N-2) -1];
   r1=roots(c1);
   r2=roots(c2);
   
   t1=sort(angle(r1)*2*sqrt(3)/pi);
   t2=sort(angle(r2)*2*sqrt(3)/pi);
   
   bsu1=t1(11);
   bsu2=t2(19);
   
   f1=abs(1-abs(sinc(N*bsu1/2/sqrt(3))./sinc(bsu1/2/sqrt(3)))-2/ASNR(i));
   f2=abs(1-abs(sinc(N*bsu2/2/sqrt(3))./sinc(bsu2/2/sqrt(3)))-2/ASNR(i));
   
   for n=1:length(u)
      bff1(i,n)=abs(bsu1/u(n));
      bff2(i,n)=abs(bsu2/u(n));
      bff(i,n)=min(bff1(i,n),bff2(i,n));
   end
end

figure
semilogy(u,bff)
xlabel('{\itu}_{\its}','Fontsize',14)
ylabel('{\itB}_{\itf}','Fontsize',14)
%title('10-element uniform linear array')
grid

FS=12;
q=length(u);
axis([0 1 0.001 100])
p=text(0.83,bff(1,q)*1.5,'ASNR=10 dB');
set(p,'Fontsize',FS)
p=text(0.83,bff(2,q)*1.5,'ASNR=20 dB');
set(p,'Fontsize',FS)
p=text(0.83,bff(3,q)*1.5,'ASNR=30 dB');
set(p,'Fontsize',FS)
p=text(0.83,bff(4,q)*1.5,'ASNR=40 dB');
set(p,'Fontsize',FS)
