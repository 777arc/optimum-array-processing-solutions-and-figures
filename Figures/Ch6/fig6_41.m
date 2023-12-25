%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.41
% Xiaomin Lu 6/26/00
% updated by K. Bell 11/30/00
% updated by Lillian Xu 02/13/2001, K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ex 6.6.7
%MPDR_DL, with array perturbation
% Optimal LNR & Array Gain vs. Array Purturbation .
%
%  ui = +/- 0.30, INR = 30dB (equal power)
%  ua = 0(Signal Matched)
  

clear all
close all

N = 10;
BWNN = 4/N;
n = (-(N-1)/2:(N-1)/2)';

LNRrange = -20:50;

SNRrange = [-5 0 10 20];

% uncomment to get figure in text
sigmaPrange = 0:0.005:0.20;
%sigmaPrange = 0:0.05:0.20;

%uncomment to get figure in text
totalTrialNumber = 500;
%totalTrialNumber = 20;

Us = 0;
Ui1 = -0.30;
Ui2 = 0.30;
INR1 = 10^(30/10);
INR2 = INR1;

C = exp(j*n*pi*[0]);								 %MPDR
f = [1];


k1 = 1;
for SNR = 10.^(SNRrange/10)
    disp(['loop ' int2str(k1) ' of ' int2str(length(SNRrange)) ' ...'])
    
    k2 = 1;
    for sigmaP = sigmaPrange
        sigmaP
        k3 = 1;
        for LNR = 10.^(LNRrange/10)
            Gain(k1,k2,k3) = 0;
            for trial = 1:totalTrialNumber
                d1 = (n/2+sigmaP*randn(N,1))*Ui1 + sigmaP*randn(N,1)*sqrt(1-Ui1^2);
                d2 = (n/2+sigmaP*randn(N,1))*Ui2 + sigmaP*randn(N,1)*sqrt(1-Ui2^2);
                ds = (n/2+sigmaP*randn(N,1))*Us  + sigmaP*randn(N,1)*sqrt(1-Us^2);
                Vi1 = exp(j*2*pi*d1);
                Vi2 = exp(j*2*pi*d2);
                Vs  = exp(j*2*pi*ds);
                Ss = SNR*Vs*Vs';
                Sn = eye(N) + INR1*Vi1*Vi1'+ INR2*Vi2*Vi2';
                Sx = Ss + Sn + LNR*eye(N);     					% Diagonal Loading
                W = inv(Sx)*C*inv(C'*inv(Sx)*C)*f;           
                SINR0 = (W'*Ss*W)/(W'*Sn*W);
                SINRi = SNR/(1+INR1+INR2);
                Gain(k1,k2,k3) = Gain(k1,k2,k3) + SINR0/SINRi;
            end	% end of trial
            k3 = k3 + 1;
        end		% end of LNR
        k2 = k2 + 1;
    end			% end of sigmaP
    k1 = k1 + 1;
end				% end of SNR

Gain = 10*log10( abs(Gain)/totalTrialNumber );

for k1 = 1:length(SNRrange)
    for k2 = 1:length(sigmaPrange)
        [A(k1,k2), I(k1,k2)] = max( Gain(k1,k2,:) );
    end
end
Load = -21 +  I;
plot(sigmaPrange, Load(1,:),'-',sigmaPrange, Load(2,:),'--');
hold on
plot(sigmaPrange, Load(3,:),'-.',sigmaPrange, Load(4,:),':');
grid
xlabel('\sigma_{\itp}/\lambda','Fontsize',14)
ylabel('Optimal {\itLNR} (dB)','Fontsize',14)
%title(['MPDR\_DL with perturbation, ui = +/-0.30(30dB each), ua = 0, ' int2str(totalTrialNumber) ' trials'])
h=legend('{\itSNR}=-5 dB','{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB',4);
set(h,'Fontsize',12)

figure
plot(sigmaPrange, A(1,:),'-',sigmaPrange, A(2,:),'--');
hold on
plot(sigmaPrange, A(3,:),'-.',sigmaPrange, A(4,:),'-*');
grid
xlabel('\sigma_{\itp}/\lambda','Fontsize',14)
ylabel('Optimal gain (dB)','Fontsize',14)
%title(['MPDR\_DL with perturbation, ui = +/-0.30(30dB each), ua = 0, ' int2str(totalTrialNumber) ' trials'])
h=legend('{\itSNR}=-5 dB','{\itSNR}=0 dB','{\itSNR}=10 dB','{\itSNR}=20 dB',4);
set(h,'Fontsize',12)
axis([0 0.2 0 50])

