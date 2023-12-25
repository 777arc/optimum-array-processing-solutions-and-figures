%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.16
% Zhi Tian 7/8/99
% Updated by K. Bell 10/13/00, 7/23/01, 10/23/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
N=10;
amf=[-(N-1)/2:(N-1)/2]';

BWnn=4/N;
du=0.001:0.1:10;

nn=kron([1:N]',ones(1,N))-kron([1:N],ones(N,1));

for i=1:length(du)
   S=sinc(nn*du(i)/2);
   A=svd(S);
   Ass(:,i)=A(1:2)/trace(S);
   
   us=du(i)/2/sqrt(3);
   Bc=sinc(N*us)/sinc(us);
   Lss(:,i)=1/2*[1+abs(Bc);1-abs(Bc)];
end
Ass=10*log10(Ass); %/max(max(Ass)));
Lss=10*log10(Lss);
figure
semilogx(du/BWnn,Ass(1,:),du/BWnn,Lss(1,:),'--',du/BWnn,Ass(2,:),du/BWnn,Lss(2,:),'--')
h=legend('Exact eigenvalues','Approximate eigenvalues',4);
set(h,'Fontsize',12)
axis([0.01,40,-50,0])
ylabel('Normalized eigenvalues (dB)','Fontsize',14)
xlabel('{\itu}_{\Delta}/{\itBW}_{\itNN}','Fontsize',14)
%title('10-element uniform linear array')
grid
