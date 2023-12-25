%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.15
% Gerry Tian, Lillian Xu 
% Updated by K. Bell 1/30/04
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

K=[1:6]*N;
u=[-1:0.01:1];
L=200;

% opt
Sn=Vi*INR*Vi'+eye(N);
Pn=Sn/(INR+1);
Rxo=Vs*SNR*Vs'+Sn;
Wo=inv(Sn)*Wqui;
Wo=Wo/sum(Wo);
Ag=abs(Wo'*Vs)^2/real(Wo'*Pn*Wo);
	
for n=1:length(K)
    vp=0;
    for l=1:L
        
        % generating matrices
        Xs=Vs*sqrt(SNR/2)*(randn(1,K(n))+j*randn(1,K(n)));
        Xi=Vi*sqrt(INR/2)*(randn(1,K(n))+j*randn(1,K(n)));
        Xw=(randn(N,K(n))+j*randn(N,K(n)))/sqrt(2);
        
        Xn=Xi+Xw;  Xx=Xn+Xs;
        Rx=Xn*Xn'/K(n);
        
        W=inv(Rx)*Wqui;
        W=W/sum(W);
        
        vp=vp+abs(W'*Vs)^2/real(W'*Pn*W);
    end
    ra(n)=10*log10(Ag/(vp/L));
end
figure
plot(K,ra,'*-.')
xlabel('Number of snapshots ({\it K})')
ylabel('Array gain reduction (dB)')
