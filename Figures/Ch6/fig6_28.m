%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.28
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 11;
n = (-(N-1)/2:(N-1)/2)';
BWNN = 4/N;

a1= -[0 0.5 0.7 0.9 0.99];
na = size(a1,2);

phi = BWNN*(0.0:0.01:2);
ni = length(phi);

INR = 10^(10/10);
Ao = zeros(na,ni);
A1 = zeros(na,ni);
for m=1:na
    for i=1:ni
        vs = exp(j*pi*n*phi(i));
        rxx0 = 1/(1-abs(a1(m))^2);  % noise power
        r = [1+abs(a1(m))^2 a1(m) zeros(1,N-2)];
        c = [1+abs(a1(m))^2 a1(m) zeros(1,N-2)];
        Sxinv = toeplitz(c,r);
        Sxinv(1,1)=1;
        Sxinv(N,N) = 1;
        Sxn = inv(Sxinv)*INR/rxx0+eye(N);
        Sxninv = inv(Sxn);
        w = inv(vs'*Sxninv*vs)*Sxninv*vs;
        A1(m,i) = (1/real(w'*inv(Sxninv)*w))*(INR+1);
    end
end
vi = ones(N,1);
Sx = INR*vi*vi'+eye(N);
Sxinv = inv(Sx);
vs = exp(j*pi*n*phi);
rho_sq = abs(vs'*vi/N).^2;
A2 = N*((1+INR)/(1+N*INR))*(1+N*INR*(1-rho_sq));

h1=plot(phi/BWNN,10*log10(A1),'-');
xlabel('{\Delta}{\itu}/{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
hold on
h2=plot(phi/BWNN,10*log10(A2),'--');

axis([0 2 0 25])

for m=1:na
    text(1.7, 10*log10(A1(m,ni-10))-1, ['{\ita}_{1}=' num2str(abs(a1(m)))],'Fontsize',12)
end

h3=plot([1 1]*BWNN/2,[0 25],'-.');
hold off
h=legend([h2 h1(1) h3],'Single plane wave','AR processes','First null',4);
set(h,'Fontsize',12)
set(gca,'XTick',[0:0.2:2])

INR = 10^(20/10);
Ao = zeros(na,ni);
A1 = zeros(na,ni);
for m=1:na
    for i=1:ni
        vs = exp(j*pi*n*phi(i));
        rxx0 = 1/(1-abs(a1(m))^2);  % noise power
        r = [1+abs(a1(m))^2 a1(m) zeros(1,N-2)];
        c = [1+abs(a1(m))^2 a1(m) zeros(1,N-2)];
        Sxinv = toeplitz(c,r);
        Sxinv(1,1)=1;
        Sxinv(N,N) = 1;
        Sxn = inv(Sxinv)*INR/rxx0+eye(N);
        Sxninv = inv(Sxn);
        w = inv(vs'*Sxninv*vs)*Sxninv*vs;
        A1(m,i) = (1/real(w'*inv(Sxninv)*w))*(INR+1);
    end
end
vi = ones(N,1);
Sx = INR*vi*vi'+eye(N);
Sxinv = inv(Sx);
vs = exp(j*pi*n*phi);
rho_sq = abs(vs'*vi/N).^2;
A2 = N*((1+INR)/(1+N*INR))*(1+N*INR*(1-rho_sq));

figure

h1=plot(phi/BWNN,10*log10(A1),'-');
xlabel('{\Delta}{\itu}/{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
hold on
h2=plot(phi/BWNN,10*log10(A2),'--');

axis([0 2 0 35])

for m=1:na
    text(1.7, 10*log10(A1(m,ni-10))-1.5, ['{\ita}_{1}=' num2str(abs(a1(m)))],'Fontsize',12)
end

h3=plot([1 1]*BWNN/2,[0 35],'-.');
hold off
h=legend([h2 h1(1) h3],'Single plane wave','AR processes','First null',4);
set(h,'Fontsize',12)
set(gca,'XTick',[0:0.2:2])