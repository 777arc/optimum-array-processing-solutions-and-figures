%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.40
% Sinc-pattern with three nulls equispaced over the sector (0.18, 0.26)
% Lillian Xiaolan Xu 3/24/99
% K. Bell 9/5/00, Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=21; 
us=0;  

u=-1:0.001:1;
n=conj(-(N-1)/2:(N-1)/2)';
vs=exp(i*n.*pi.*us); 

nulls=[0.21 0.22 0.23];         

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

b=0;     
for m=1:N
    b=b+Wo(m)*exp(-i*(-(N+1)/2+m)*pi*u);
end;
bdb1=20*log10(abs(b));
b1=real(b);


figure
plot(u,bdb1);
hold on
axis([-1 1 -80 10])
line([nulls(2),nulls(2)],[-80 10]);
grid on
xlabel('\it u','Fontsize',16)
ylabel('Beam pattern (dB)','Fontsize',16)
