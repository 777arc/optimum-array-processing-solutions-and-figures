%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.37
% K. Bell 7/26/01, 11/13/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ex 6.6.5

clear all
close all

N = 10;
SNR = 10.^([0:10:30]/10);

sigmag = 10.^([-70:1:0]/10);

for s=1:4
    A = SNR(s)*(N+sigmag)./((1+(N-1)*sigmag*SNR(s)).^2 + (N-1)*(N+sigmag).*sigmag*SNR(s)^2);
    switch s
    case 1
        p = '-';
    case 2
        p='--';
    case 3
        p='-.';
    case 4
        p=':';
    end
    plot(10*log10(sigmag),10*log10(A),p)
    hold on
    
end
h=legend('{\itSNR}_{\itin}=0 dB','{\itSNR}_{\itin}=10 dB','{\itSNR}_{\itin}=20 dB','{\itSNR}_{\itin}=30 dB',3);
set(h,'Fontsize',12)
xlabel('Variance of gain error (dB)','Fontsize',14)
ylabel('Average {\itSNR}_{\ito} (dB)','Fontsize',14)
