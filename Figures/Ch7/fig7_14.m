%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.14
% Gerry Tian, Lillian Xu 
% Updated by K. Bell 1/30/04
% Note: not same as figure in text
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
nu = length(u);
V=exp(j*amf*pi*u);

L=500;
for n=1:length(K)
    vp=0;
    for l=1:L
        
        % generating matrices
        Xs=Vs*sqrt(SNR/2)*(randn(1,K(n))+j*randn(1,K(n)));
        Xi=Vi*sqrt(INR/2)*(randn(1,K(n))+j*randn(1,K(n)));
        Xw=(randn(N,K(n))+j*randn(N,K(n)))/sqrt(2);
        
        Xn=Xi+Xw;  
        Xx=Xn+Xs;
        Rx=Xn*Xn'/K(n);
        
        W=inv(Rx)*Wqui;
        W=W/sum(W);
        bp=abs(W'*V).^2;
        
        peaks = [find(diff([sign(diff(bp(nu-1:nu))) sign(diff(bp))])==-2)] ;  % local maxima on grid, wrap in u
        lp = length(peaks);
        [mbu,mbp]=min(abs(u(peaks))); % closest to MRA
        I=find([1:1:lp]~=mbp);
        peaks=peaks(I);               % remove main beam peak
        [bpeak,I]=sort(bp(peaks));
        SL = bpeak(lp-1);             % maximum of remaining peaks  
        vp=vp+10*log10(SL);           % average in dB
    end
    sll(n)=vp/L;
end
figure
plot(K,sll,'*-.')
xlabel('Number of snapshots ({\it K})')
ylabel('Highest sidelobe (dB)')
