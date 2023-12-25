% ex7.3.8 no signal mismatch
% fig 7.21w-z

N=10;

To=2/N;

amf=[-(N-1)/2:(N-1)/2]';
ua=0.03;
Va=exp(j*pi*amf*ua);
Vs=ones(N,1);
SNR=10.^(10/10);

ui=[0.29 0.45];D=2;
INR=10.^(20/10);
Vi=exp(j*amf*pi*ui);
Sn=INR*Vi*Vi'+eye(N);
Pn=Sn/(1+2*INR);

Mc=1;
C=Vs;
g=Vs'*C/N;

%Mc=3;
%C=[Vs j*amf.*Vs -amf.^2.*Vs];
%g=Vs'*C/N;

B=orth(eye(N)-C*inv(C'*C)*C');
Wq=C*inv(C'*C)*g';

BWnn=4/N;
alfa=real(To-Wq'*Wq);
Ks=[1:10 20 50 100 200 500 1000]*N;
L=200;

for i=1:length(Ks)
   K=Ks(i);
   ARvl=0; ARrl=0; Betal=0; Tvl=0; Trl=0;
   for l=1:L
		% generating data matrices
		Xs=Va*sqrt(SNR/2)*(randn(1,K)+j*randn(1,K));
		Xi=Vi*sqrt(INR/2)*(randn(D,K)+j*randn(D,K));
		Xw=(randn(N,K)+j*randn(N,K))/sqrt(2);
      Xn=Xi+Xw;  Xx=Xn+Xs;
		Rx=Xx*Xx'/K;
       
      Rz=B'*Rx*B;
      Wa=inv(Rz)*B'*Rx*Wq;	% czo: approx VL
      War=Wa;
   	cq=real(Wa'*Wa)-alfa;
	   lmta=0;
      if cq>0
         V_czo=inv(Rz)*Wa;
         bq=2*real(Wa'*V_czo);
         aq=real(V_czo'*V_czo);
         dq=bq^2-4*aq*cq;    
         lmta=(bq-sqrt(max(0,dq)))/2/aq;
         Wa=Wa-lmta*V_czo;
         War=inv(Rz+lmta*eye(N-Mc))*B'*Rx*Wq;
      end 
      Wv=Wq-B*Wa;
      %Wv=gsc_vl(N,K,Xx,C,B,g',To);
      %[Wv,lmta]=direct(N,Rx,To,C,g');
      ARvl=ARvl+real(abs(Wv'*Va)^2/(Wv'*Pn*Wv));
      Betal=Betal+abs(lmta);
      Tvl=Tvl+Wv'*Wv;
      
      Wr=Wq-B*War;	% czo: approx VL
      ARrl=ARrl+real(abs(Wr'*Va)^2/(Wr'*Pn*Wr));
      Trl=Trl+Wr'*Wr;
   
   end
   Wvbp(:,i)=Wv;
   Wrbp(:,i)=Wr;
   ARv(i)=10*log10(ARvl/L);
   ARr(i)=10*log10(ARrl/L);
   Beta(i)=10*log10(Betal/L);
   Tv(i)=Tvl/L;
   Tr(i)=Trl/L;
end

u=[-1:0.001:1];
Bp10=20*log10(Wvbp(:,1)'*exp(j*pi*amf*u));
Bp100=20*log10(Wvbp(:,10)'*exp(j*pi*amf*u));
Bp1000=20*log10(Wvbp(:,13)'*exp(j*pi*amf*u));
Bpr10=20*log10(Wrbp(:,1)'*exp(j*pi*amf*u));
Bpr100=20*log10(Wrbp(:,10)'*exp(j*pi*amf*u));
Bpr1000=20*log10(Wrbp(:,13)'*exp(j*pi*amf*u));

SNR=10*log10(SNR);
INR=10*log10(INR);

save ex7_3_8 L Tv Tr ARv ARr Beta SNR INR ua ui Ks To u Bp10 Bp100 Bp1000 Bpr10 Bpr100 Bpr1000 Wrbp Wvbp

figure
semilogx(Ks,ARv,Ks,ARr,'--')
legend('VL-SMI-1','VL-SMI-2')
xlabel('Number of Snapshots: K')
ylabel('Array Gain (dB)')
title(['Array Gain vs. K: SNR=10dB, INR=20dB, ui=[0.29 0.45], ua=',num2str(ua),', To=',num2str(To)]) 
grid

figure
semilogx(Ks, Beta)
xlabel('Number of Snapshots: K')
ylabel('Amount of Variable Loading: \beta (dB)')
title(['Variable loading vs. K: SNR=10dB, INR=20dB, ui=[0.29 0.45], ua=',num2str(ua),', To=',num2str(To)]) 
grid

figure
subplot(311)
plot(u,Bp10,u,Bpr10,'--')
legend('VL-SMI-1','VL-SMI-2')
hold on
plot([ui(1),ui(1)],[-40,10],'-.')
plot([ui(2),ui(2)],[-40,10],'-.')
plot([ua,ua],[-40,10],'--')
axis([-1,1,-40,10])
xlabel('u')
ylabel('K = N')
title(['Representative Beampattern for K=10: SNR=10dB, INR=20dB, To=',num2str(To)]) 
grid

subplot(312)
plot(u,Bp100,u,Bpr100,'--')
legend('VL-SMI-1','VL-SMI-2')
hold on
plot([ui(1),ui(1)],[-40,10],'-.')
plot([ui(2),ui(2)],[-40,10],'-.')
plot([ua,ua],[-40,10],'--')
axis([-1,1,-40,10])
xlabel('u')
ylabel('K = 10 N')
%title(['Representative Beampattern for K=100: SNR=10dB, INR=20dB, To=',num2str(To)]) 
grid

subplot(313)
plot(u,Bp1000,u,Bpr1000,'--')
legend('VL-SMI-1','VL-SMI-2')
hold on
plot([ui(1),ui(1)],[-40,10],'-.')
plot([ui(2),ui(2)],[-40,10],'-.')
plot([ua,ua],[-40,10],'--')
axis([-1,1,-50,10])
xlabel('u')
ylabel('K = 100 N')
%title(['Representative Beampattern for K=1000: SNR=10dB, INR=20dB, To=',num2str(To)]) 
grid