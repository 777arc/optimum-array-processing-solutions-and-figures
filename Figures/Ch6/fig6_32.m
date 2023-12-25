%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.32
% MVDR with Signal Mismatch
% Array Gain Ratio 
% Lillian Xiaolan Xu 04/09/2001
% Updated by K. Bell 7/26/01, 11/12/01
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
for ui = [0.3]
    for INR = 10.^([10]/10)
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
    plot(ua,R)
    xlabel('{\itu}_{\ita}','Fontsize',14)
    ylabel('{\itA}_{\itmvdr} /{\itA}_{\ito}({\bfv}_{\itm})','Fontsize',14)
    %title(['N=',num2str(N),', u_{I}=',num2str(ui)])
    axis([-0.1 0.1 0.4 1])
    grid
    
    k0=k0+1;
end

