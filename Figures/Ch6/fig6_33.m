%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.33
% MVDR with Signal Mismatch
% Array Gain Ratio 
% Lillian Xiaolan Xu 09/21/2000, 04/09/2001
% updated by K. Bell 7/26/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
Vm = ones(N,1);
u = -1:1/1000:1;
ua=-0.1:0.001:0.1;
SNR = 10.^([30]/10);
k0 = 1;
k1 = 1;
Sso=SNR*Vm*Vm';
figure
for ui = [0.0433 0.02]
    for INR = 10.^([10 20 30]/10)
        inr(k1)=INR;
        Vi = exp(j*n*pi*ui);
        Sn = INR*Vi*Vi' + eye(N);
        W = inv(Sn)*Vm/(Vm'*inv(Sn)*Vm);
        for num=1:size(ua,2)
            Va = exp(j*n*pi*ua(num));
            Ss = SNR*Va*Va';
            A=W'*Ss*W/(W'*Sn*W)/(SNR/(INR+1));
            Ao=W'*Sso*W/(W'*Sn*W)/(SNR/(INR+1));
            R(num,k1)=real(A/Ao);
        end
        k1 = k1 + 1;
    end
    subplot(2,1,k0)
    plot(ua,R(:,3*(k0-1)+1),ua,R(:,3*(k0-1)+2),'--',ua,R(:,3*(k0-1)+3),'-.')
    ylabel('{\itA}_{\itmvdr} /{\itA}_{\ito}({\bfv}_{\itm})    ','Fontsize',12)
    %title(['N=',num2str(N),', u_{I}=',num2str(ui)])
    grid
    
    k0=k0+1;
end

subplot(2,1,1)
xlabel('(a)','Fontsize',12)
h=legend(['{\itINR}=',num2str(10*log10(inr(1))),' dB'],...
    ['{\itINR}=',num2str(10*log10(inr(2))),' dB'],...
    ['{\itINR}=',num2str(10*log10(inr(3))),' dB']);
set(h,'Fontsize',12)
axis([-0.1 0.1 0 5])

subplot(2,1,2)
xlabel('{\itu}_{\ita}','Fontsize',12)
h=legend(['{\itINR}=',num2str(10*log10(inr(1))),' dB'],...
    ['{\itINR}=',num2str(10*log10(inr(2))),' dB'],...
    ['{\itINR}=',num2str(10*log10(inr(3))),' dB']);
set(h,'Fontsize',12)
axis([-0.1 0.1 0 20])
text(-0.003,-5.5,'(b)','Fontsize',12)
