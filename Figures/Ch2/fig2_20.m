%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.20
% Effect of element spacing on beam pattern
% Xin Zhang
% Lillian Xu updated 09/2000
% updated by K. Bell 7/22/2001, 10/4/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2).';
psi = pi*(-3:0.001:3);
w = 1/N*[ones(N,1)];
d = [1/4 1/2 1];

figure
for k =1:length(d)
    vv = exp(j*2*d(k)*n*psi);
    B(k,:) = 20*log10(abs(w'*vv));
    subplot(3,1,k)
    plot(psi/pi,B(k,:));
    hold on
    plot([-1 -1],[-25 0],'--');
    plot([1 1],[-25 0],'--');
    axis([-3 3 -25 5])
    set(gca,'YTick',[-25 -20 -15 -10 -5 0])
    grid on
    xlabel('\itu')
    if k == 1
        plot([-1 1],[2.5 2.5])
        plot([-1 -0.9],[2.5 0.5])
        plot([-1 -0.9],[2.5 4.5])
        plot([0.9 1],[0.5 2.5])
        plot([0.9 1],[4.5 2.5])
        text(-0.4,7,'Visible region')
    end
    hold off
end
