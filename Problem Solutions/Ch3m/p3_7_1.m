%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.7.3
% Initial 40dB Chebychev pattern with two 1st & 2nd order nulls at 0.235,0.265.
% K. Bell 9/24/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

N=41;
D = 0.5*[-(N-1)/2:1:(N-1)/2].';
u=[-1:0.001:1];
v = exp(j*2*pi*D*u);

sidelobe=40;   	% *dB below the main lobe maximum
R=10^(sidelobe/20);
x0=cosh(1/(N-1)*acosh(R));
% chebychev weights
Wd = poly(exp(j*2*acos(cos((2*[1:1:N-1]-1)*pi/(2*(N-1)))/x0))).';

Wd=Wd./sum(Wd);

% Example 3.7.3 null constraints
nulls=[0.22 0.24 0.26 0.28];             
C = exp(j*2*pi*D*nulls);
Pc=C*inv(C'*C)*C';
Wo=(eye(N)-Pc)*Wd;

Bo = Wo'*v;

% 0,1,and 2 order nulls
nulls=[0.236 0.265];             
vn = exp(j*2*pi*D*nulls);
dvn = (j*2*pi*D*ones(1,length(nulls))).*vn;
dv2n = ((j*2*pi*D*ones(1,length(nulls))).^2).*vn;
C4 = [vn dvn];
Pc=C4*inv(C4'*C4)*C4';
W1=(eye(N)-Pc)*Wd;
B1 = W1'*v;

C6 = [vn dvn dv2n];
Pc=C6*inv(C6'*C6)*C6';
W2=(eye(N)-Pc)*Wd;
B2 = W2'*v;

subplot(2,1,1)
h1=plot(u,20*log10(abs(Bo)));
hold on
h2=plot(u,20*log10(abs(B1)),'--r');
hold off
axis([-1 1 -80 10]);
legend([h1 h2],'0-order nulls (4 const.)','0 and 1-order nulls (4 const.)',1)
xlabel('\it u')
ylabel('Beam pattern (dB)')
title('Problem 3.7.1')
grid on

subplot(2,1,2)

h1=plot(u,20*log10(abs(Bo)));
hold on
h3=plot(u,20*log10(abs(B2)),'--r');
hold off
axis([-1 1 -80 10]);
legend([h1 h3],'0-order nulls (4 const.)','0,1,and 2-order nulls (6 const.)',1)
xlabel('\it u')
ylabel('Beam pattern (dB)')
grid on
