%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7.8, 7.9, 7.10
% Gerry Tian 6/24/99
% Updated by K. Bell 1/30/04
% Note: Figure 7.9 not same as text (text actually used LNR=15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N=10;
SNR=10.^([0:10:30]/10);
INR=10^(20/10);

K=[0.2 0.4 0.8 1 2 3 4 6 10 20 50 100]*N;
L=500;
J=rot90(eye(N));

uo=0; ui=[0.29 0.45]; D=2;
amf=[-(N-1)/2:(N-1)/2]';
Vs=exp(j*pi*amf*uo);
Vi=exp(j*pi*amf*ui);
LNR=10.^(20/10);

for i=1:length(SNR)
    
    % SINR,optimal in steady state: MVDR, MPDR 
    Sn=Vi*INR*Vi'+eye(N);
    Rxo=Vs*SNR(i)*Vs'+Sn;
    
    Wpo=inv(Vs'*inv(Rxo)*Vs)*Vs'*inv(Rxo);
    SINRpo=SNR(i)*abs(Wpo*Vs)^2/real(Wpo*Sn*Wpo');
    SINRpopt(i)=10*log10(abs(SINRpo));
    
    Wvo=inv(Vs'*inv(Sn)*Vs)*Vs'*inv(Sn);
    SINRvo=SNR(i)*abs(Wvo*Vs)^2/real(Wvo*Sn*Wvo');
    SINRvopt(i)=10*log10(abs(SINRvo));
    
    for n=1:length(K)
        
        % SMI: MPDR, MVDR
        Sp=0;Spb=0;Sv=0;Svb=0;
        for l=1:L
            % generating data matrices
            Xs=Vs*sqrt(SNR(i)/2)*(randn(1,K(n))+j*randn(1,K(n)));
            Xi=Vi*sqrt(INR/2)*(randn(D,K(n))+j*randn(D,K(n)));
            Xw=(randn(N,K(n))+j*randn(N,K(n)))/sqrt(2);
            
            Xn=Xi+Xw;  Xx=Xn+Xs;
            
            Rx=Xx*Xx'/K(n)+LNR*eye(N);
            Rxfb=(Rx+J*conj(Rx)*J)/2;
            
            Rn=Xn*Xn'/K(n)+LNR*eye(N);
            Rnfb=(Rn+J*conj(Rn)*J)/2;
            
            % Weighting
            Wmpdr=inv(Vs'*inv(Rx)*Vs)*Vs'*inv(Rx);
            Wmpfb=inv(Vs'*inv(Rxfb)*Vs)*Vs'*inv(Rxfb);
            
            Wmvdr=inv(Vs'*inv(Rn)*Vs)*Vs'*inv(Rn);
            Wmvfb=inv(Vs'*inv(Rnfb)*Vs)*Vs'*inv(Rnfb);
            
            % SINR,smi
            Sp =Sp +SNR(i)*abs(Wmpdr*Vs)^2/real(Wmpdr*Sn*Wmpdr');
            Spb=Spb+SNR(i)*abs(Wmpfb*Vs)^2/real(Wmpfb*Sn*Wmpfb');
            Sv =Sv +SNR(i)*abs(Wmvdr*Vs)^2/real(Wmvdr*Sn*Wmvdr');
            Svb=Svb+SNR(i)*abs(Wmvfb*Vs)^2/real(Wmvfb*Sn*Wmvfb');
            
        end
        
        SINRp(i,n) =10*log10(abs(Sp/L));
        SINRpb(i,n)=10*log10(abs(Spb/L));
        SINRv(i,n) =10*log10(abs(Sv/L));
        SINRvb(i,n)=10*log10(abs(Svb/L));
    end
end

SNR=10*log10(SNR);
INR=10*log10(INR);
LNR=10*log10(LNR);

%%% figure 7.8 %%%

figure
semilogx(K,SINRp(1,:)-SINRpopt(1),K,SINRp(2,:)-SINRpopt(2),'-.',K,SINRp(3,:)-SINRpopt(3),'--',K,SINRp(4,:)-SINRpopt(4),':');
hold on
semilogx(K,SINRpb(1,:)-SINRpopt(1),'x-',K,SINRpb(2,:)-SINRpopt(2),'o-.',K,SINRpb(3,:)-SINRpopt(3),'+--',K,SINRpb(4,:)-SINRpopt(4),'*:');
legend('{\it SNR} = 0 dB','{\it SNR} = 10 dB','{\it SNR} = 20 dB','{\it SNR} = 30 dB','{\it SNR} = 0 dB, FB','{\it SNR} = 10 dB, FB','{\it SNR} = 20 dB, FB','{\it SNR} = 30 dB, FB',4)
grid on
xlabel('Number of snapshots (K)')
ylabel('Normalized average {\it SINR_0} (dB)')
axis([10 1000 -50 0])


%%% figure 7.10 %%%
figure
semilogx(K,SINRp(1,:)-SINRpopt(1),K,SINRp(2,:)-SINRpopt(2),'-.',K,SINRp(3,:)-SINRpopt(3),'--',K,SINRp(4,:)-SINRpopt(4),':');
hold on
semilogx(K,SINRpb(1,:)-SINRpopt(1),'x-',K,SINRpb(2,:)-SINRpopt(2),'o-.',K,SINRpb(3,:)-SINRpopt(3),'+--',K,SINRpb(4,:)-SINRpopt(4),'*:');
axis([1,1000,-50,0])
grid
legend('{\it SNR} = 0 dB','{\it SNR} = 10 dB','{\it SNR} = 20 dB','{\it SNR} = 30 dB','{\it SNR} = 0 dB, FB','{\it SNR} = 10 dB, FB','{\it SNR} = 20 dB, FB','{\it SNR} = 30 dB, FB',4)
grid on
xlabel('Number of snapshots (K)')
ylabel('Normalized average {\it SINR_0} (dB)')
axis([2 1000 -50 0])


%%% figure 7.9 %%%
figure
semilogx(K,SINRv(1,:)-SINRvopt(1),K,SINRvb(1,:)-SINRvopt(1),'-.');
legend('MVDR','MVDR-FB',4)
xlabel('Number of Snapshots (K)')
ylabel('Normalized Average SINR (dB)')
title(['MVDR.SMI beamformers: INR = ',num2str(INR),' dB, \sigma^{2}_{L} = 20 dB'])
grid
axis([10,1000,-5,0])