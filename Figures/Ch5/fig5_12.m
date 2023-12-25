%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.12
% K. Bell 10/13/00, 7/23/01, 10/23/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigenvalues for 2 uncorrelated sources

clear all
close all

N=20;
a=10.^([0 3 10 20]/10);      % S2/S1
BWNN = 4/N;
del_u = 10.^([-2:0.01:0.8]);    % del_u/BWNN
Bc12 = sinc(N*BWNN*del_u/2)./sinc(BWNN*del_u/2);
figure
for n = 1:4
   l1=0.5*(1+a(n)+sqrt((1-a(n))^2+4*a(n)*abs(Bc12).^2));
   l2=0.5*(1+a(n)-sqrt((1-a(n))^2+4*a(n)*abs(Bc12).^2));
   semilogx(del_u,10*log10([l1;l2]))
   hold on
end
hold off
axis([0.01,3,-35,25])
xlabel('{\it\Deltau}/{\itBW}_{\itNN}','Fontsize',14)
ylabel('Normalized eigenvalues (dB)','Fontsize',14)
grid
