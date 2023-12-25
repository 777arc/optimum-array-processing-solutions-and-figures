%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6.27
%  K. Bell 7/25/01, 11/12/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 11;
n = (-(N-1)/2:(N-1)/2)';
BWNN = 4/N;

a1= -[0.2 0.5 0.7];
na = size(a1,2);

phi = BWNN*(0.01:0.01:2);
ni = length(phi);

Ao = zeros(na,ni);
A1 = zeros(na,ni);
for m=1:na
    for i=1:ni
        vs = exp(j*pi*n*phi(i));
        r = [1+abs(a1(m))^2 a1(m) zeros(1,N-2)];
        c = [1+abs(a1(m))^2 a1(m) zeros(1,N-2)];
        Sxinv = toeplitz(c,r);
        Sxinv(1,1)=1;
        Sxinv(N,N) = 1;
        w = inv(vs'*Sxinv*vs)*Sxinv*vs;
        rxx0 = 1/(1-abs(a1(m))^2);
        A1(m,i) = (1/real(w'*inv(Sxinv)*w))*rxx0;
        Ao(m,i) = N*(1+2*((N-1)/N)*(abs(a1(m))/(1-abs(a1(m))^2))*(abs(a1(m))-cos(pi*phi(i))));
    end
end

plot(phi/BWNN,10*log10(Ao),'-');
hold on
plot(phi/BWNN,10*log10(A1),'--')
xlabel('{\Delta}{\itu}/{\itBW}_{\itNN}','Fontsize',14)
ylabel('Array gain (dB)','Fontsize',14)
axis([0 2 0 20])

for m=1:na
    text(1.7, 10*log10(Ao(m,ni)), ['{\ita}_{1}=' num2str(abs(a1(m)))],'Fontsize',12)
end

hold off