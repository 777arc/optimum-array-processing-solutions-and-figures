%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.11
% Gerry Tian 6/24/99
% Updated by K. Bell 1/30/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=10;
SNR=10.^(10/10);
INR=10^(20/10);
LNRdb=[-10:10:20];
LNR=10.^(LNRdb/10);

K=4;
L=1;
J=rot90(eye(N));

uo=0; ui=[0.29 0.45]; D=2;
amf=[-(N-1)/2:(N-1)/2]';
Vs=exp(j*pi*amf*uo);
Vi=exp(j*pi*amf*ui);

u=[-1:0.01:1];
V=exp(j*pi*amf*u);
% steady state
for i=1:length(LNR)
    
    Xs=Vs*sqrt(SNR/2)*(randn(1,K)+j*randn(1,K));
    Xi=Vi*sqrt(INR/2)*(randn(D,K)+j*randn(D,K));
    Xw=(randn(N,K)+j*randn(N,K))/sqrt(2);
    
    Xn=Xi+Xw;  Xx=Xn+Xs;
    
    Rx=Xx*Xx'/K+LNR(i)*eye(N);
    Rxfb=(Rx+J*conj(Rx)*J)/2;
    
    Wp=inv(Rxfb)*Vs*inv(Vs'*inv(Rxfb)*Vs);
    Bp(i,:) =20*log10(abs(Wp'*V));
end

save ex7_3_4e Bp LNR LNRdb N K ui u 
figure
for i=1:length(LNR)
    %figure
    subplot(2,2,i)
    plot(u,Bp(i,:))
    hold on
    plot([ui(1),ui(1)],[-70,10],'--')
    plot([ui(2),ui(2)],[-70,10],'--')
    xlabel('u')
    ylabel('Beam pattern (dB)')
    title(['{\it N} = 10,{\it K} = 4, \sigma_{\itL}^{2} = ',num2str(LNRdb(i)),' dB'])
    grid
    if i==1
        axis([-1,1,-50,20])
    else
        axis([-1,1,-60,10])
    end
end

