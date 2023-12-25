%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.37
% Illustration of the Parks-McClellan algorithm
% for equiripple approximation
% Lillian Xu 3/24/99
% K. Bell 9/22/00, K. Bell 7/22/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

up=0.3;
us=0.4;

u1=[0.001:0.002:up];
u2=[up+0.001:0.002:1.001];
u=[u1 u2];
L=7;
dp=0.1;
ds=0.1;

psi_set=[0.06 0.1 0.18 up us 0.5 0.6 0.8 1]*pi;

bd=[1 1 1 1 0 0 0 0 0];
w=[[1 1 1 1]*ds/dp [1 1 1 1 1]];
bd_s=[ones(1,size(u1,2)) zeros(1,size(u2,2))];
w_s=[ones(1,size(u1,2))*ds/dp ones(1,size(u2,2))];
error=w_s*bd_s';
%    while error>0.1

x=cos(psi_set);

for k=1:L+2
    b(k)=1;
    for i=1:L+2
        if i~=k
            b(k)=b(k)/(x(k)-x(i));
        end
    end
end

d1=0;
d2=0;
for k=1:L+2
    d1=d1+b(k)*bd(k);
    d2=d2+b(k)*(-1)^(k+1)/w(k);
end

delta=d1/d2;


for k=1:L+1
    f(k)=bd(k)-(-1)^(k+1)*delta/w(k);
    d(k)=b(k)/(x(k)-x(L+2));
end

b1=0;
b2=0;
for k=1:L+1
    b1=b1+d(k)*f(k)./(cos(u*pi)-x(k));
    b2=b2+d(k)./(cos(u*pi)-x(k));
end

bo=b1./b2;

nb = length(bo);
k1 = [find(diff([0 sign(diff(bo))])==-2)];   % local maxima on grid
k2 = [find(diff([0 sign(diff(bo))])==2)];   % local minima on grid
if ((sign(diff(bo(1:2))))==-1)     % decreasing at startpoint
    k1=[1 k1];
else
    k2 = [1 k2];
end
if ((sign(diff(bo(nb-1:nb))))==1)  % increasing at endpoint
    k1=[k1 nb];
else
    k2 = [k2 nb];
end

new_set=[u(k1) u(k2) up us];
%error=w_s*(bd_s-bo)';
% end

plot([up us 1],[1+delta -delta bo(k2(4))],'o')
hold on
plot(new_set,[bo(k1) bo(k2) 1+delta -delta],'*');
plot(u,bo);
plot([0 1],[0 0],'--');
plot([0 1],[1 1],'--');
plot([0 psi_set(4)/pi],[1-delta 1-delta],'--');
plot([psi_set(5)/pi 1],[-delta -delta],'--');
plot([0 psi_set(4)/pi],[1+delta 1+delta],'--');
plot([psi_set(5)/pi 1],[delta delta],'--');

psi_set=[0.06 0.1 0.18 0.5 0.6 0.8]*pi;
k=1;
for num=1:size(u,2)
    if k<L
        if abs(psi_set(k)/pi-u(num))<=0.0015
            plot(psi_set(k)/pi,bo(num),'o');
            k=k+1;
        end
    end
end

plot([up up],[0 1],':')
plot([us us],[0 1],':')
h=legend('Initial set of {\it u(i)}','New set of {\it u(i)}');
hold off
axis([0 1 -0.2 1.2])
text(up-0.02,-0.04,'\it u_p','Fontsize',12)
text(us-0.02,-0.07,'\it u_s','Fontsize',12)
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)
