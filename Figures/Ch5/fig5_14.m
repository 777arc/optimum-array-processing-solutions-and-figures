%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.14
% Zhi Tian 7/8/99
% Updated by K. Bell 10/12/00, 7/23/01, 10/23/01
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
   Ass(:,i)=A/trace(S);
end
Ass=10*log10(Ass); %/max(max(Ass)));
figure
semilogx(du/BWnn,Ass)
%semilogx(du/BWnn,Ass(1,:),du/BWnn,Lss(1,:),'--',du/BWnn,Ass(2,:),du/BWnn,Lss(2,:),'--')
axis([0.01,40,-50,0])
ylabel('Normalized eigenvalues (dB)','Fontsize',14)
xlabel('{\itu}_{\Delta}/{\itBW}_{\itNN}','Fontsize',14)
%title('10-element uniform linear array')
grid
