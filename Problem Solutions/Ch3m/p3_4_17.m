%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.4.17
% K. Bell 9/29/00
% Updated by K. Bell 9/02/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
N=[10:10:40];
s = ['-k';'--';':k';'-.'];

d=0.5;
SL = [-20:-1:-100];

ns = length(SL);
BWNN = zeros(4,ns);
for m = 1:4
   for q=1:ns
      R = 10^(-SL(q)/20);
      x0 = cosh(acosh(R)/(N(m)-1));
      BWNN(m,q) = (2/(pi*d))*acos(cos(pi/(2*(N(m)-1)))/x0);
   end
   h(m)=plot(SL,BWNN(m,:),s(m,:));
   hold on

end
hold off
legend(h,['N=' int2str(N(1))],['N=' int2str(N(2))],['N=' int2str(N(3))],['N=' int2str(N(4))])
xlabel('Sidelobe Level (dB)')
ylabel('BWNN (u)')
title('Problem 3.4.17')