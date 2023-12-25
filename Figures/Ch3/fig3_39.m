%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.39 null steering
% (a) Initial sinc-pattern
% (b) Pattern with a null of zero order imposed at u=0.22
% (c) With null of first order too
% (d) With null of second order too
% 
% Lillian Xiaolan Xu 3/24/99
% Last updated by K. Bell 7/22/01, 9/30/01, 10/17/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=21; 
us=0;     
u=-1:0.001:1;
n=conj(-(N-1)/2:(N-1)/2)';
vs=exp(i*n.*pi.*us); 

% ----------------null constraints
nulls=[0.22];   
vnull=exp(i*n.*pi.*nulls(1)); 
C=vnull;
if size(nulls,2)>1
   for m=2:size(nulls,2)
      C=[C exp(i*n.*pi.*nulls(m))]; 
   end;
end;

C1=[C i*n.*vnull];
C2=[C1 -n.*n.*vnull];

Pc=C*inv(C'*C)*C';
Pc1=C1*inv(C1'*C1)*C1';
Pc2=C2*inv(C2'*C2)*C2';

Wd=1/N*ones(N,1);
Wo=(Wd'*(eye(N)-Pc))';
Wo1=(Wd'*(eye(N)-Pc1))';
Wo2=(Wd'*(eye(N)-Pc2))';
b=0;     
for m=1:N
   b=b+Wd(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb1=20*log10(abs(b));
b1=real(b);
b=0;     
for m=1:N
   b=b+Wo(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb2=20*log10(abs(b));
b2=real(b);
b=0;     
for m=1:N
   b=b+Wo1(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb3=20*log10(abs(b));
b3=real(b);
b=0;     
for m=1:N
   b=b+Wo2(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb4=20*log10(abs(b));
b4=real(b);


figure
subplot(2,2,1)
plot(u,bdb1);
hold on
axis([-1 1 -80 10])
title('(a) Uniform','Fontsize',14)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
line([nulls(1),nulls(1)],[-80 10]);
grid on

subplot(2,2,2)
plot(u,bdb2);
hold on
axis([-1 1 -80 10])
title('(b) Zero-order null','Fontsize',14)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
line([nulls(1),nulls(1)],[-80 10]);
grid on

subplot(2,2,3)
plot(u,bdb3);
hold on
axis([-1 1 -80 10])
title('(c) First-order null','Fontsize',14)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
line([nulls(1),nulls(1)],[-80 10]);
grid on

subplot(2,2,4)
plot(u,bdb4);
hold on
axis([-1 1 -80 10])
title('(d) Second-order null','Fontsize',14)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
line([nulls(1),nulls(1)],[-80 10]);
grid on
