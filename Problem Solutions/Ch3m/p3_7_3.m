%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.7.3
% K. Bell 9/23/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

N1 = 11;
N2 = 5;
d=0.5;

D = [d*[-(N2-1)/2:1:(N2-1)/2]-7.5 d*[-(N1-1)/2:1:(N1-1)/2] d*[-(N2-1)/2:1:(N2-1)/2]+7.5].';
N=length(D);

subplot(3,1,1)
h1=plot(D.', zeros(1,N),'x');


N3 = 35;
D3 = d*[-(N3-1)/2:1:(N3-1)/2].';
hold on
h2=plot(D3.', zeros(1,N3),'o');
xlabel('z/\lambda')
title('Element positions')
legend([h1 h2],'21-element array','35-element SLA')
hold off

u = [-1:0.001:1];
v = exp(j*2*pi*D*u);
v3 = exp(j*2*pi*D3*u);

w_unf = ones(N,1)/N;
Bd = w_unf'*v;

B3 = ones(N3,1)'*v3/N3;

subplot(3,1,2)
h1=plot(u,20*log10(abs(Bd)));
hold on
h2=plot(u,20*log10(abs(B3)),'--');
axis([-1 1 -40 0])
grid on
xlabel('u')
ylabel('Beam pattern (dB)')
legend([h1 h2],'21-element array','35-element SLA')
title('Problem 3.7.3(a)')


u_null = 0.5/N;
C = exp(j*2*pi*D*u_null);
Pc_orth = eye(N)-C*inv(C'*C)*C';
wc = Pc_orth*w_unf;

Bc = wc'*v;
subplot(3,1,3)
h1=plot(u,20*log10(abs(Bd)),'--');
hold on
h2=plot(u,20*log10(abs(Bc)),'-');
plot(u_null*[1 1],[-60 0],'g')
hold off
axis([-1 1 -40 0])
grid on
xlabel('u')
ylabel('Beam pattern (dB)')
legend([h1 h2],'Desired pattern','Null Constraint')
title('Problem 3.7.3(b)')
