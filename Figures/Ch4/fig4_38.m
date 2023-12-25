%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.38
% Beampatterns for a uniform circular array of 20 elements with uniform weighting:
%     phi=0 deg and phi=90 deg
% Lillian Xiaolan Xu 9/20/99
% Updated by K. Bell 9/30/00
% Updated by Lillian 12/06/2000, K. Bell 7/23/01
% Function called: polardb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
N=20; 

lim=-40;
uy_set=-1:0.005:1;
ny = size(uy_set,2);
ux_set=0;
nx = size(ux_set,2);
b = zeros(nx,ny);

for y=1:ny
	uy=uy_set(y);
	for x=1:nx
		ux=ux_set(x);

		p=sqrt(ux*ux+uy*uy);
		e=acos(ux/p);
		
		b(x,y)=0;
		for m=-N:N
			if m<0
			b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*(-1)^(m*N)*besselj(abs(m*N),10*p);
			else
			b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
			end
		end
	end
end
b=b/max(max(b));


bdb=20*log10(abs(b));

for x=1:nx
	for y=1:ny
		if bdb(x,y)<lim
			bdb(x,y)=lim;
		end
		if sqrt(ux_set(x)^2+uy_set(y)^2)>1
			bdb(x,y)=lim;
		end;
	end
end

bdb1=bdb;


figure
polardb(asin(uy_set),bdb1,lim)
hold on
polardb(asin(uy_set)+pi,bdb1,lim)
hold off
