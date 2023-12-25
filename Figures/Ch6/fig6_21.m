%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.21
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
vs = ones(N,1);

u = -1:0.01:1;
v = exp(j*n*pi*u);
SNR = 1;
INR = 10.^([20]/10);

ui = [0.3 0.18 0.02];

for i=1:3
    VI = [exp(j*n*pi*ui(i)) exp(-j*n*pi*ui(i))];
    R = INR*VI*VI'+eye(N);
    
    w = inv(vs'*inv(R)*vs)*inv(R)*vs;
    B = w'*v;
    
    figure
    plot(u,10*log10(abs(B).^2),'-')
    hold on
    if i==3
        plot(ui(i)*[1 1],[-10 20],'--')
        plot(-ui(i)*[1 1],[-10 20],'--')
        axis([-1 1 -10 30])
    else
        plot(ui(i)*[1 1],[-30 -10],'--')
        plot(-ui(i)*[1 1],[-30 -10],'--')
        axis([-1 1 -30 0])
    end
    xlabel('{\itu}','Fontsize',14)
    ylabel('Beam pattern (dB)','Fontsize',14)
    grid on
end