%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.23
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 21;
n = (-(N-1)/2:(N-1)/2)';
vs = ones(N,1);

u = -1:0.001:1;
v = exp(j*n*pi*u);
SNR = 1;
INR = 10.^([0 10 20]/10);

ui = [0.21 0.22 0.23];

for i=1:3
    VI = [exp(j*n*pi*ui)];
    R = INR(i)*VI*VI'+eye(N);
    
    w = inv(vs'*inv(R)*vs)*inv(R)*vs;
    B = w'*v;
    
    figure
    plot(u,10*log10(abs(B).^2),'-')
    hold on
    plot(ui(1)*[1 1],[-50 -10],'--')
    plot(ui(2)*[1 1],[-50 -10],'--')
    plot(ui(3)*[1 1],[-50 -10],'--')
    axis([-1 1 -100 0])
    xlabel('{\itu}','FontSize',14)
    ylabel('Beam pattern (dB)','FontSize',14)
    grid on
end