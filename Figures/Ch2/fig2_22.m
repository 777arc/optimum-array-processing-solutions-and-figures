%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2.22
% Effect of element spacing on beam pattern
% Lillian Xu
% Last updated by K. Bell 7/22/01, 10/4/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all

N=10;
n=(-(N-1)/2:(N-1)/2)';
w=ones(N,1)/N;
u=-3:0.01:3;

d=[2/3 1/2];
theta=[30 90];
figure
s = ['(a)';'(b)'];
for m=1:2
    beam=(w.*exp(i*2*pi*n*d(m)*sin(theta(m)/180*pi)))'...
        *exp(i*2*pi*n*d(m)*u);
    subplot(2,1,m);
    plot(u,20*log10(abs(beam)));
    hold on
    axis([-3 3 -25 5])
    ylabel('Beam pattern (dB)')
    xlabel('\itu')
    text(0,-32,s(m,:),'HorizontalAlignment','center')
%title(['N=',num2str(N),', d=',num2str(d(m),2),' lambda, theta-Bar=',num2str(theta(m)),' degrees'])
    grid on
end

