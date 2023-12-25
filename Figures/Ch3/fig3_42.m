%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.42
% Initial 40dB Chebychev pattern with eight nulls equispaced over the sector(0.22,0.36). Sidelobe cancellation=51dB
% Lillian Xiaolan Xu 3/24/99
% Last updated by K. Bell 9/22/00, 7/22/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

N=41;
u=-1:0.001:1;
n=-(N-1)/2:(N-1)/2;

sidelobe=40;   	% *dB below the main lobe maximum
R=10^(sidelobe/20);
x0=cosh(1/(N-1)*acosh(R));
% chebychev weights
Wd = poly(exp(j*2*acos(cos((2*[1:1:N-1]-1)*pi/(2*(N-1)))/x0))).';
Wd=Wd./sum(Wd);

% ----------------null constraints

nulls=[0.22:0.02:0.36];     
n=conj(-(N-1)/2:(N-1)/2)';
vnull=exp(i*n.*pi.*nulls(1)); 
C=vnull;
if size(nulls,2)>1
   for m=2:size(nulls,2)
      C=[C exp(i*n.*pi.*nulls(m))]; 
   end;
end;
Pc=C*inv(C'*C)*C';
Wo=(Wd'*(eye(N)-Pc))';

b=0;     
for m=1:N
   b=b+Wo(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb1=20*log10(abs(b));
b=0;     
for m=1:N
   b=b+Wd(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb2=20*log10(abs(b));

figure

plot(u,bdb1);
axis([-1 1 -80 10]);
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
grid on
hold on

hold off
return
