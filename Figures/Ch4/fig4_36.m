%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.36
% Beampatterns for a uniform circular array:
%     (a) vertical theta0=0, phi=0,pi
%     (b) vertical theta0=0, phi=pi/2,3*pi/2
% Lillian Xiaolan Xu 9/20/99	
% last updated 9/30/00 by K. Bell, 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           --------------------------                                    %
%                          Circular Ring array beampatterns                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%% (a)
N=7; 
lim=-60;
ux_set=-1:0.005:1;
uy_set=0;

for x=1:size(ux_set,2)
	ux=ux_set(x);
	for y=1:size(uy_set,2)
		uy=uy_set(y);
		b0(x,y)=besselj(0,10*ux);
		p=sqrt(ux*ux+uy*uy);
		e=acos(ux/p);		
		b(x,y)=0;
		for m=-N:N
			b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
		end
	end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
	for y=1:size(uy_set,2)
		if bdb(x,y)<lim
			bdb(x,y)=lim;
		end
		if sqrt(ux_set(x)^2+uy_set(y)^2)>1
			bdb(x,y)=lim;
		end;
	end
end

bdb1=bdb;
%%%%%%%%%%%%%%%%%%%

N=10;
for x=1:size(ux_set,2)
	ux=ux_set(x);
	for y=1:size(uy_set,2)
		uy=uy_set(y);
		p=sqrt(ux*ux+uy*uy);
		e=acos(ux/p);
		b(x,y)=0;
		for m=-N:N
			b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
		end
	end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
	for y=1:size(uy_set,2)
		if bdb(x,y)<lim
			bdb(x,y)=lim;
		end
		if sqrt(ux_set(x)^2+uy_set(y)^2)>1
			bdb(x,y)=lim;
		end;
	end
end

bdb2=bdb;
%%%%%%%%%%%%%%%%%%%

N=20;
for x=1:size(ux_set,2)
	ux=ux_set(x);
	for y=1:size(uy_set,2)
		uy=uy_set(y);
		p=sqrt(ux*ux+uy*uy);
		e=acos(ux/p);
		b(x,y)=0;
		for m=-N:N
			b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
		end
	end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
	for y=1:size(uy_set,2)
		if bdb(x,y)<lim
			bdb(x,y)=lim;
		end
		if sqrt(ux_set(x)^2+uy_set(y)^2)>1
			bdb(x,y)=lim;
		end;
	end
end
bdb3=bdb;
%%%%%%%%%%%%%%%%%%%

N=25;
for x=1:size(ux_set,2)
	ux=ux_set(x);
	for y=1:size(uy_set,2)
		uy=uy_set(y);
		p=sqrt(ux*ux+uy*uy);
		e=acos(ux/p);
		b(x,y)=0;
		for m=-N:N
			b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
		end
	end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
	for y=1:size(uy_set,2)
		if bdb(x,y)<lim
			bdb(x,y)=lim;
		end
		if sqrt(ux_set(x)^2+uy_set(y)^2)>1
			bdb(x,y)=lim;
		end;
	end
end

figure
subplot(2,1,1)
plot(ux_set,bdb1,'-',ux_set,bdb2,'--',ux_set,bdb3,'-.',ux_set,bdb,':')
axis([-1 1 lim 0])
h=legend('{\it N}=7','{\it N}=10','{\it N}=20','{\it N}=25');
set(h,'Fontsize',12)
xlabel(['\it u_x'],'Fontsize',14); %  ( theta=',num2str(theta),' degrees)'])
ylabel('Beam pattern (dB)','Fontsize',14)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (b)

clear
N=7;
lim=-60;
uy_set=-1:0.005:1;
ux_set=0;

for y=1:size(uy_set,2)
    uy=uy_set(y);
    for x=1:size(ux_set,2)
        ux=ux_set(x);
        b0(x,y)=besselj(0,10*uy);
        p=sqrt(ux*ux+uy*uy);
        e=acos(ux/p);
        b(x,y)=0;
        for m=-N:N
            b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
        end
    end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
    for y=1:size(uy_set,2)
        if bdb(x,y)<lim
            bdb(x,y)=lim;
        end
        if sqrt(ux_set(x)^2+uy_set(y)^2)>1
            bdb(x,y)=lim;
        end;
    end
end

bdb1=bdb;
%%%%%%%%%%%%%%%%%%%

N=10;
for x=1:size(ux_set,2)
    ux=ux_set(x);
    for y=1:size(uy_set,2)
        uy=uy_set(y);
        p=sqrt(ux*ux+uy*uy);
        e=acos(ux/p);
        b(x,y)=0;
        for m=-N:N
            b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
        end
    end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
    for y=1:size(uy_set,2)
        if bdb(x,y)<lim
            bdb(x,y)=lim;
        end
        if sqrt(ux_set(x)^2+uy_set(y)^2)>1
            bdb(x,y)=lim;
        end;
    end
end
bdb2=bdb;
%%%%%%%%%%%%%%%%%%%

N=20;
for x=1:size(ux_set,2)
    ux=ux_set(x);
    for y=1:size(uy_set,2)
        uy=uy_set(y);
        p=sqrt(ux*ux+uy*uy);
        e=acos(ux/p);
        b(x,y)=0;
        for m=-N:N
            b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
        end
    end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
    for y=1:size(uy_set,2)
        if bdb(x,y)<lim
            bdb(x,y)=lim;
        end
        if sqrt(ux_set(x)^2+uy_set(y)^2)>1
            bdb(x,y)=lim;
        end;
    end
end
bdb3=bdb;
%%%%%%%%%%%%%%%%%%%

N=25;
for x=1:size(ux_set,2)
    ux=ux_set(x);
    for y=1:size(uy_set,2)
        uy=uy_set(y);
        p=sqrt(ux*ux+uy*uy);
        e=acos(ux/p);        
        b(x,y)=0;
        for m=-N:N
            b(x,y)=b(x,y)+j^(m*N)*exp(-i*m*N*e)*besselj(m*N,10*p);
        end
    end
end
b=b/max(max(b));
bdb=20*log10(abs(b));

for x=1:size(ux_set,2)
    for y=1:size(uy_set,2)
        if bdb(x,y)<lim
            bdb(x,y)=lim;
        end
        if sqrt(ux_set(x)^2+uy_set(y)^2)>1
            bdb(x,y)=lim;
        end;
    end
end

subplot(2,1,2)
plot(uy_set,bdb1,uy_set,bdb2,'--',uy_set,bdb3,'-.',uy_set,bdb,':')
axis([-1 1 lim 0])
h=legend('{\it N}=7','{\it N}=10','{\it N}=20','{\it N}=25');
set(h,'Fontsize',12)
xlabel(['\it u_y'],'Fontsize',14); %  ( theta=',num2str(theta),' degrees)'])
ylabel('Beam pattern (dB)','Fontsize',14)
grid on
