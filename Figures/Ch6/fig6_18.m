%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.18
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 10;
n = (-(N-1)/2:(N-1)/2)';
BWNN = 4/N;
vs = ones(N,1);

ui = BWNN*(0:0.01:2.75);
ni = length(ui);

SNR = 1;

INR = 10.^([20]/10);
vi = exp(j*n*pi*ui);
rho_sq = abs(vs'*vi/N).^2;
Ao = N*((1+INR)/(1+N*INR))*(ones(1,ni)+N*INR*(ones(1,ni)-rho_sq));
SINR = SNR/(INR+1);

sidelobe=20;   	% *dB below the main lobe maximum
R=10^(sidelobe/20);
x0=cosh(1/(N-1)*acosh(R));
% chebychev weights
wdc = poly(exp(j*2*acos(cos((2*[1:1:N-1]-1)*pi/(2*(N-1)))/x0))).';
wdc = wdc/sum(wdc);

ADC = (SNR/SINR)./(INR*abs(wdc'*vi).^2 + real(wdc'*wdc));
h1=plot(ui/BWNN,10*log10(ADC),'-');
hold on
h2=plot(ui/BWNN,10*log10(Ao./ADC),'--');

plot([0.23 0.23], [0 35],':')
plot([0.5 0.5], [0 35],':')
xlabel('{\itu}_{\itI} /{\itBW}_{\itNN}','FontSize',14)
ylabel('dB','FontSize',14)

h=legend([h1 h2],'{\itA}_{\itDC},20 dB ripple','{\itA}_{\ito} /{\itA}_{\itDC}',0);
set(h,'FontSize',12)
axis([0 3 0 35])
hold off