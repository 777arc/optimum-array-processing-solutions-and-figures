%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.10
% Zhi Tian 7/7/99
% Updated by K. Bell 10/12/00, 7/23/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=10;
amf=[-(N-1)/2:(N-1)/2]';

BWnn=4/N;
du=0.001:0.01:2;

nn=kron([1:N]',ones(1,N))-kron([1:N],ones(N,1));

for i=1:length(du)
   us=du(i);
   Bc=sin(N*pi*us/2)/N/sin(pi*us/2);
   Lss(:,i)=[1+abs(Bc);1-abs(Bc)];
end
Lss=10*log10(Lss);
figure
semilogx(du/BWnn,Lss)
axis([0.01,3,-35,5])
xlabel('{\it\Deltau}/{\itBW}_{\itNN}','Fontsize',14)
ylabel('Normalized eigenvalues (dB)','Fontsize',14)
grid
