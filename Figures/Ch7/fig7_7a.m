%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.7(a)
% Gerry Tian 11/29/98
% Updated by K. Bell 1/30/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=10;
SNR=10.^([0:10:30]/10);
INR=10^(20/10);

K=[1 2 3 4 6 10 20 50 100]*N;
L=500;
J=rot90(eye(N));

uo=0; ui=[0.29 0.45]; D=2;
amf=[-(N-1)/2:(N-1)/2]';
Vs=exp(j*pi*amf*uo);
Vi=exp(j*pi*amf*ui);

%LCMV
C=[Vs (j*amf).*Vs -(amf.^2).*Vs];
g=C'*ones(N,1)/N;
Mc=3;
Pvt=eye(N)-C*inv(C'*C)*C';
B=orth(Pvt);

for i=1:length(SNR)
    
    % SINR,optimal in steady state: MVDR, MPDR 
    Sn=Vi*INR*Vi'+eye(N);
    Rxo=Vs*SNR(i)*Vs'+Sn;
    
    Wvo=inv(Vs'*inv(Sn)*Vs)*Vs'*inv(Sn);
    SINRvo=SNR(i)*abs(Wvo*Vs)^2/(Wvo*Sn*Wvo');
    SINRvopt(i)=10*log10(abs(SINRvo));
    
    tmp=inv(Sn)*C;
    Wgsco=tmp*inv(C'*tmp)*g;
    
    SINRgo=SNR(i)*abs(Wgsco'*Vs)^2/real(Wgsco'*Sn*Wgsco);
    SINRgopt(i)=10*log10(abs(SINRgo));
    
    for n=1:length(K)
        %Analytic: mvdr
        Ev(i,n)=10*log10((K(n)+2-(N-2))/(K(n)+1)*SINRgo);
        
        % SMI: LCMV
        Sv=0;
        for l=1:L
            % generating data matrices
            Xs=Vs*sqrt(SNR(i)/2)*(randn(1,K(n))+j*randn(1,K(n)));
            Xi=Vi*sqrt(INR/2)*(randn(D,K(n))+j*randn(D,K(n)));
            Xw=(randn(N,K(n))+j*randn(N,K(n)))/sqrt(2);
            
            Xn=Xi+Xw;  Xx=Xn+Xs;
            
            Rn=Xn*Xn'/K(n);
            
            % Weighting
            tmp=inv(Rn)*C;
            Wgsc=tmp*inv(C'*tmp)*g;
            
            % SINR,smi
            Sv =Sv + SNR(i)*abs(Wgsc'*Vs)^2/real(Wgsc'*Sn*Wgsc);
        end
        
        SINRgsc(i,n) =10*log10(abs(Sv/L));
        
    end
end

SNR=10*log10(SNR);
INR=10*log10(INR);

figure
semilogx(K,SINRgsc(1,:)-SINRgopt(1),K,SINRgsc(2,:)-SINRgopt(2),'-.',K,SINRgsc(3,:)-SINRgopt(3),'--',K,SINRgsc(4,:)-SINRgopt(4),':');
hold on
semilogx(K,Ev(1,:)-SINRgopt(1),'x',K,Ev(2,:)-SINRgopt(2),'x',K,Ev(3,:)-SINRgopt(3),'x',K,Ev(4,:)-SINRgopt(4),'x');
axis([10 1000 -5 0])
legend('{\it SNR} = 0 dB','{\it SNR} = 10 dB','{\it SNR} = 20 dB','{\it SNR} = 30 dB','{\it E}[\rho]',4)
xlabel('Number of snapshots ({\it K})')
ylabel('Normalized average {\it SINR_0} (dB)')
grid on
