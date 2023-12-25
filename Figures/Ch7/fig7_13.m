%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.13
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

L=200;
K=(1:6)*N;

for n=1:length(K)
	sp=0;
	for l=1:L
		Xi=Vi*sqrt(INR/2)*(randn(1,K(n))+j*randn(1,K(n)));
		Xw=(randn(N,K(n))+j*randn(N,K(n)))/sqrt(2);
		Xn=Xi+Xw;  
		Rn=Xn*Xn'/K(n);
		Rw=Xw*Xw'/K(n);
		lmta=svd(Rn);
		sp=sp+(lmta(2)/lmta(10))^2;
	end
	spd(n)=10*log10(sp/L);
end

figure
plot(K,spd,'*-')
axis([10,60,0,100])
xlabel('Number of snapshots ({\it K})')
ylabel('Eigenvalue Spread (dB)')
