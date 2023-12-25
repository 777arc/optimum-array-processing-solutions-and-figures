%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.5.5 (a), (c), (e)
% K. Bell 10/6/99
% Updated 9/02/02 by K. Bell
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

d=0.5;
u = [-1:0.01:1];
uu = [0:0.0001:1];
psi_0 = 1/sqrt(2)
N=21;
D = d*[-(N-1)/2:1:(N-1)/2].';
a = psi_0*sinc([0:1:(N-1)/2]*psi_0);
Bu = a(1)+2*a(2:(N+1)/2)*cos([1:1:(N-1)/2].'*u*pi);
BBu = a(1)+2*a(2:(N+1)/2)*cos([1:1:(N-1)/2].'*uu*pi);
del1 = max(BBu);
del2 = abs(min(BBu));
del(1) = max(del1-1,del2);
[p,i] = find(BBu<=1-del(1));
psi_p = uu(min(i));
[p,i] = find(BBu>=del(1));
psi_s = uu(max(i));
del_psi(1) = psi_s-psi_p;

figure(1)
subplot(3,1,1)
h1=plot(u,Bu,'k');
hold on
plot([-1 -psi_0 -psi_0 psi_0 psi_0 1],[0 0 1 1 0 0])
wu = [flipud(a.');a(2:(N+1)/2).'];

Whamm = 0.54+0.46*cos([0:1:(N-1)/2]*2*pi/N);
B = a(1)*Whamm(1)+2*(a(2:(N+1)/2).*Whamm(2:(N+1)/2))*cos([1:1:(N-1)/2].'*u*pi);
BB = a(1)*Whamm(1)+2*(a(2:(N+1)/2).*Whamm(2:(N+1)/2))*cos([1:1:(N-1)/2].'*uu*pi);
del1 = max(BB);
del2 = abs(min(BB));
del(2) = max(del1-1,del2);
[p,i] = find(BB<=1-del(2));
psi_p = uu(min(i));
[p,i] = find(BB>=del(2));
psi_s = uu(max(i));
del_psi(2) = psi_s-psi_p;
wh = [flipud((a.').*(Whamm(1:(N+1)/2).')); (a(2:(N+1)/2).').*(Whamm(2:(N+1)/2).')];

A = -20*log10(del(2));
if A<21
   beta = 0;
elseif A>=21 & A<=50
   beta = 0.5842*(A-21)^(0.4)+0.07886*(A-21);
else
   beta = 0.1102*(A-8.7);
end
beta
Wk = besseli(0,beta*sqrt(1-(2*[0:1:(N-1)/2]/(N-1)).^2));
Wk = Wk/Wk(1);   
Bk = a(1)+2*(a(2:(N+1)/2).*Wk(2:(N+1)/2))*cos([1:1:(N-1)/2].'*u*pi);
BB = a(1)+2*(a(2:(N+1)/2).*Wk(2:(N+1)/2))*cos([1:1:(N-1)/2].'*uu*pi);
del1 = max(BB);
del2 = abs(min(BB));
del(3) = max(del1-1,del2);
[p,i] = find(BB<=1-del(3));
psi_p = uu(min(i));
[p,i] = find(BB>=del(3));
psi_s = uu(max(i));
del_psi(3)= psi_s-psi_p;
wk = [flipud((a.').*(Wk(1:(N+1)/2)).'); (a(2:(N+1)/2).').*(Wk(2:(N+1)/2).')];

h2=plot(u,B,'--g');
h3=plot(u,Bk,'-.r');
hold off
legend([h1 h2 h3],'uniform','Hamming','Kaiser')
%xlabel('u=\psi/\pi')
ylabel('Beampattern')
title(['Problem 3.5.5(a,c,e), Hamming, N = ' int2str(N)])
axis([-1 1 -0.1 1.1])
hold off

subplot(3,1,2)
h1=plot(u,Bu,'k');
hold on
plot([-1 -psi_0 -psi_0 psi_0 psi_0 1],[0 0 1 1 0 0])
h2=plot(u,B,'--g');
h3=plot(u,Bk,'-.r');
hold off
legend([h1 h2 h3],'uniform','Hamming','Kaiser')
%xlabel('u=\psi/\pi')
ylabel('Beampattern')
hold off
axis([-1 1 0.95 1.05])
subplot(3,1,3)
h1=plot(u,Bu,'k');
hold on
plot([-1 -psi_0 -psi_0 psi_0 psi_0 1],[0 0 1 1 0 0])
h2=plot(u,B,'--g');
h3=plot(u,Bk,'-.r');
hold off
legend([h1 h2 h3],'uniform','Hamming','Kaiser')
xlabel('u=\psi/\pi')
ylabel('Beampattern')
hold off
axis([-1 1 -0.05 0.05])

set(gcf,'Paperposition',[0.25 1 8 9])

w = [wu wh wk]
del
del_psi
