%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.15
% Zhi Tian 7/7/99
% Updated by K. Bell 10/12/00, 7/23/01, 10/23/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=10;
amf=[-(N-1)/2:(N-1)/2]';

BWnn=4/N;
du=0.2;

nn=kron([1:N]',ones(1,N))-kron([1:N],ones(N,1));

S=sinc(nn*du/2);
[U,A,Ut]=svd(S);
W=U(:,1:2);
W=W/sum(W(:,1));

u=[-1:0.01:1];
Bp=20*log10(abs(W'*exp(j*amf*pi*u)));
figure
plot(u,Bp(1,:),u,Bp(2,:),'--')
axis([-1,1,-50,0])
xlabel('\itu','Fontsize',14)
ylabel('Eigenbeams (dB)','Fontsize',14)
%title([num2str(N),'-element ULA, u_{\Delta}=',num2str(du)])
grid
