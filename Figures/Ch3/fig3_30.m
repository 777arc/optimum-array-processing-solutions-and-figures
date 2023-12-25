%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.30
% Villeneuve roots
% Lillian Xu
% Modified by Xin Zhang 3/24/99
% Updated by K. Bell 9/5/00
% Last updated by Lillian Xiaolan Xu 09/20/2000, K. Bell 7/23/01, 9/30/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
N=21;
M=(N-1)/2;

n1=6;
d=1/2;

sidelobe=20;
R=10.^(sidelobe/20);
x0=cosh(1/2/M*acosh(R));
u=-1:0.001:1;

lim=-80;

%-------------------Aperture

v=N*d*u;
theta=acos(u);
a=acosh(R)/pi;
aa=a.*a;
m1=1;
m2=1;
for num=1:n1-1
	vn=n1.*sqrt((aa+(num-0.5)^2)./(aa+(n1-0.5)^2));
	m1=m1.*(1-v.*v./vn./vn);
	m2=m2.*(1-v.*v/num/num);
end;
b=sinc(v).*m1./m2;

for num=1:size(b,2)
	if b(num)==Inf | b(num)==-Inf
		if num~=1 & num~=size(b,2)
			b(num)=(b(num+1)+b(num-1))/2;
		elseif num==1
			b(num)=b(num+1);
		else
			b(num)=b(num-1);
		end
	end
end

bdb=20*log10(abs(b));

%----------------   Array   --------

b1=sinc(N*u*d)./sinc(u*d);

for num=1:n1-1
	un=1/N/d*n1.*sqrt((aa+(num-0.5)^2)./(aa+(n1-0.5)^2));
	b1=b1.*(u-un).*(u+un)./(u-num/N/d)./(u+num/N/d); 
end

for num=1:size(b1,2)
	if b1(num)==Inf | b1(num)==-Inf
		if num~=1 & num~=size(b,2)
			b1(num)=(b1(num+1)+b1(num-1))/2;
		elseif num==1
			b1(num)=b1(num+1);
		else
			b1(num)=b1(num-1);
		end
	end
end
b1=b1./max(b1);

%--------------------
c=pi*n1/N/acos(1/x0*cos((2*n1-1)*pi/4/M));

b2=sinc(N/2*u)./sinc(u/2);
for num=1:n1-1
	up=2*acos(1/x0*cos((2*num-1)*pi/4/M))/pi;
	b2=b2.*sin(pi/2*(u-c*up)).*sin(pi/2*(u+c*up))./sin(pi/2*(u-2*num/N))./sin(pi/2*(u+2*num/N));
end
for num=1:size(b2,2)
	if b2(num)==Inf | b2(num)==-Inf
		if num~=1 & num~=size(b2,2)
			b2(num)=(b2(num+1)+b2(num-1))/2;
		elseif num==1
			b2(num)=b2(num+1);
		else
			b2(num)=b2(num-1);
		end
	end
end

b2=b2./max(b2);

%------------------------
num=1:M;
uroot=2*num/N;
croot=2/pi*acos(1/x0*cos((2*num-1)*pi/4/M));
vroot=[c*croot(1:n1-1) uroot(n1:M)];

%wu=ones(N,1)/N;
%zu=roots(wu);
%wc=Dcheb(N,sidelobe);
%zc=roots(wc);
zu=exp(i*pi*uroot);
zc=exp(i*pi*croot);
zv=exp(i*pi*vroot);

figure

theta=2*pi*[0:0.001:1];

p1=1.5;
p2=2;
p3=2.5;

plot(cos(theta),sin(theta),'--');
hold on
plot(p1*cos(theta),p1*sin(theta),'--'); 
plot(p2*cos(theta),p2*sin(theta),'--'); 

plot(real(zu),imag(zu),'o');
plot(p1*real(zv),p1*imag(zv),'o');
plot(p2*real(zc),p2*imag(zc),'o');
for k=1:size(zv,2)
   plot([0 p2*1.1*real(zv(k))],[0 p2*1.1*imag(zv(k))]);
end
plot([-p2*1.2 p2*1.2],[0 0]); 
plot([0 0],[-p2*1.2 p2*1.2]); 
plot([-0.02 0],[p2*1.1 p2*1.2]); 
plot([0.02 0],[p2*1.1 p2*1.2]); 
plot([p2*1.1 p2*1.2],[-0.02 0]); 
plot([p2*1.1 p2*1.2],[0.02 0]); 


axis([-p3 p3 -p3 p3])
axis('square')
axis off
grid on
hold off
title('Roots in {\it z}-plane, {\it N}=21, {\it n}=6','Fontsize',14)
text(-0.5,-p2-0.1,'Chebychev roots','Fontsize',12)
text(-0.5,-p1-0.1,'Villeneuve roots','Fontsize',12)
text(-0.5,-1.1,'Uniform roots','Fontsize',12)
text(p2+0.2,-0.15,'Real','Fontsize',14)
text(0.1,p2+0.3,'Imaginary','Fontsize',14)
text(-0.1,-0.1,'0')
text(1.05,2.76,'-')