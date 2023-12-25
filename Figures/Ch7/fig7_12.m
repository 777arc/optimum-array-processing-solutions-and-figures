%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.12
% K. Bell 1/30/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=10;
amf=[-(N-1)/2:(N-1)/2]';

us=0;
Vs=exp(j*pi*amf*us);
SNR=10^(20/10);
ui=0.29;
Vi=exp(j*pi*amf*ui);
INR=10^(30/10);

K=([2 6])*N;

% Dolph-Chebychev weights
R = 10^(30/20);
x0 = cosh(acosh(R)/(N-1));
z = exp(j*2*acos(cos((2*[1:1:N-1]-1)*pi/(2*(N-1)))/x0));
wt = poly(z).';
Wdq = wt/sum(wt);

C=Vs;
f=1;
Pct=eye(N)-C*inv(C'*C)*C';
Wqui=Pct*Wdq + C*inv(C'*C)*f;

u=[-1:0.01:1];
nu = length(u);
V=exp(j*amf*pi*u);

for n=1:length(K)
    Xs=Vs*sqrt(SNR/2)*(randn(1,K(n))+j*randn(1,K(n)));
    Xi=Vi*sqrt(INR/2)*(randn(1,K(n))+j*randn(1,K(n)));
    Xw=(randn(N,K(n))+j*randn(N,K(n)))/sqrt(2);
    
    Xn=Xi+Xw;  
    Xx=Xn+Xs;
    Rx=Xn*Xn'/K(n);
    W=inv(Rx)*Wqui;
    W=W/sum(W);
    bp(n,:)=10*log10(abs(W'*V).^2);
end

figure
plot(u,bp(1,:),'-')
hold on
plot(u,bp(2,:),'--')
axis([-1 1 -60 10])
hold off
legend('{\it K}=2{\itN}','{\it K}=6{\itN}')
xlabel('{\it u}')
ylabel('Beam pattern (dB)')
grid on