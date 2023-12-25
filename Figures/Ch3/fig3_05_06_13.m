%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.5
% DPSS psi0=0.1pi,0.2pi, and 0.4pi
% Figure 3.6
% Beampattern of a standard 11-element linear array as a function
% Figure 3.13
% Zero plots for DPSS(0.1pi,0.2pi,0.4pi) weightings
%
% Lillian Xiaolan Xu
% Updated 4/14/99
% Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
% Functions called: DPSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
N=11; 

us=0;                           % should change us between 0 and 0.2
u=-1:0.001:1;
n=conj(-(N-1)/2:(N-1)/2)';
vs=exp(i*n.*pi.*us); 

%------------------------ DPSS Weighting
psi0=[0.1 0.2 0.4]*pi;
w1=DPSS(psi0(1),N);     		% 
w2=DPSS(psi0(2),N);			%
w3=DPSS(psi0(3),N);

z1=roots(w1);
z2=roots(w2);
z3=roots(w3);

w1=w1/max(w1);
w2=w2/max(w2);
w3=w3/max(w3);

m=-5:5;
figure
h1 = plot(m,w1);
hold on;
h2 = plot(m,w1,'o');
h3 = plot(m,w2,'--');
h4 = plot(m,w2,'o');
h5 = plot(m,w3,'-.');
h6 = plot(m,w3,'o');
hold off
h=legend([h1 h3 h5],'{\it \psi}_0=0.1\pi','{\it \psi}_0=0.2\pi','{\it \psi}_0=0.4\pi');
set(h,'Fontsize',12)
%title('DPSS for 11 elements');
xlabel('Sensors','Fontsize',14);
ylabel('Weighting','Fontsize',14);
grid on

%-----------------------------------------------------------%
% ---------------- Zero Plot.
figure
%%%%%%%%%%%%      plot 1

subplot(1,2,1);
plot(real(z1),imag(z1),'o');
axis([-1.2 1.2 -1.2 1.2])
axis('square')
grid on;
title('{\it \psi}_0=0.1\pi','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
hold off

%%%%%%%%%%%%	plot 2


subplot(1,2,2);
plot(real(z3),imag(z3),'o');
grid on;
axis([-1.2 1.2 -1.2 1.2])
axis('square')
title('{\it \psi}_0=0.4\pi','Fontsize',14);
xlabel('Real','Fontsize',14)
ylabel('Imaginary','Fontsize',14)
theta=2*pi*[0:0.01:1];
hold on
plot(cos(theta),sin(theta),'--'); 
hold off

w1=w1./sum(w1);
w2=w2./sum(w2);
w3=w3./sum(w3);

b=0;     
for m=1:N
    b=b+w1(m)*exp(i*(-(N+1)/2+m)*pi*u);
end;
bdb1=20*log10(abs(b));
b1=real(b); 
b=0;     
for m=1:N
    b=b+w2(m)*exp(i*(-(N+1)/2+m)*pi*u);
end;
bdb2=20*log10(abs(b));
b2=real(b);
b=0;     
for m=1:N
    b=b+w3(m)*exp(i*(-(N+1)/2+m)*pi*u);
end;
bdb3=20*log10(abs(b));
b3=real(b);

figure
plot(u,bdb1,u,bdb2,'--',u,bdb3,'-.');
hold on
h=legend('{\it \psi}_0=0.1\pi','{\it \psi}_0=0.2\pi','{\it \psi}_0=0.4\pi');
axis([-1 1 -80 0])
grid on
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
