%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.10.1 
% K. Bell 10/10/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
N=32;
d=0.5;
D = d*[-(N-1)/2:1:(N-1)/2].';
u = [-1:0.001:1];
v = exp(j*2*pi*D*u);
nu = length(u);
u2 = [-3:1:3]*2/N;
n2 = length(u2);
v2 = exp(j*2*pi*D*u2);

Bbs = v2/sqrt(N);
vbs = Bbs'*v;

wd = Bbs'*exp(j*2*pi*D*0)/N;
BB = wd'*vbs;

%(a)
C = Bbs'*exp(j*2*pi*D*3/N);
Pc = eye(n2)-C*inv(C'*C)*C';
wn = Pc*wd;
BN = wn'*vbs;

figure(1)
h1=plot(u, 10*log10(abs(BB.^2)),':');
hold on
h2=plot(u, 10*log10(abs(BN.^2)),'-');
h3 = plot([1 1]*3/N,[-60 0]);
hold off
axis([-1 1 -60 0])
legend([h1 h2],'B_d','B_d w/ null at u=3/32',2)
title('Problem 3.10.1(a)')
ylabel('Beampattern (dB)')
xlabel('u')

%(b)
C = Bbs'*exp(j*2*pi*D*7/N);
Pc = eye(n2)-C*inv(C'*C)*C';
wn = Pc*wd;
BN = wn'*vbs;

figure(2)
subplot(2,1,1)
h1=plot(u, 10*log10(abs(BB.^2)),':');
hold on
h2=plot(u, 10*log10(abs(BN.^2)),'-');
h3 = plot([1 1]*7/N,[-60 0]);
hold off
axis([-1 1 -60 0])
legend([h1 h2],'B_d','B_d w/ null at u=7/32',2)
title('Problem 3.10.1(b)')
ylabel('Beampattern (dB)')
xlabel('u')

C = Bbs'*exp(j*2*pi*D*13/N);
Pc = eye(n2)-C*inv(C'*C)*C';
wn = Pc*wd;
BN = wn'*vbs;

subplot(2,1,2)
h1=plot(u, 10*log10(abs(BB.^2)),':');
hold on
h2=plot(u, 10*log10(abs(BN.^2)),'-');
h3 = plot([1 1]*13/N,[-60 0]);
hold off
axis([-1 1 -60 0])
legend([h1 h2],'B_d','B_d w/ null at u=13/32',2)
ylabel('Beampattern (dB)')
xlabel('u')

%(c)
wd = Bbs'*exp(j*2*pi*D*1/N)/N;
BB = wd'*vbs;

C = Bbs'*exp(j*2*pi*D*5/N);
Pc = eye(n2)-C*inv(C'*C)*C';
wn = Pc*wd;
BN = wn'*vbs;

figure(3)
h1=plot(u, 10*log10(abs(BB.^2)),':');
hold on
h2=plot(u, 10*log10(abs(BN.^2)),'-');
h3 = plot([1 1]*3/N,[-60 0]);
hold off
axis([-1 1 -60 0])
legend([h1 h2],'B_d','B_d w/ null at u=5/32',2)
title('Problem 3.10.1(c)')
ylabel('Beampattern (dB)')
xlabel('u')