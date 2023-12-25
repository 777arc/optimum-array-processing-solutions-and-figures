%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.22
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
BWNN = 4/N;
vs = ones(N,1);

ui = BWNN*(0.01:0.01:2.5);
ni = length(ui);

SNR = 1;
INR = 10.^([-20 -10 0 10 20]/10);
nn = length(INR);

Ao = zeros(nn,ni);
for m=1:nn
    for i=1:ni
        VI = [exp(j*n*pi*ui(i)) exp(-j*n*pi*ui(i))];
        R = INR(m)*VI*VI'+eye(N);
        SINR = SNR/(2*INR(m)+1);
        
        w = inv(vs'*inv(R)*vs)*inv(R)*vs;
        Ao(m,i) = SNR/(N*SINR*real(w'*R*w));
    end
end

plot(ui/BWNN,10*log10(Ao),'-');
xlabel('{\itu}_{\itI} /{\itBW}_{\itNN}','FontSize',14)
ylabel('{\itA}_{\ito} /{\itN} (dB)','FontSize',14)
axis([0 2.5 -10 25])

for m=1:nn
    if m==1
        text(1.5, 10*log10(Ao(m,ni))-1, ['{\itINR}=' num2str(10*log10(INR(m))) ' dB'],'FontSize',12)
    else
        text(1.5, 10*log10(Ao(m,ni))+1, ['{\itINR}=' num2str(10*log10(INR(m))) ' dB'],'FontSize',12)
    end
end

hold off