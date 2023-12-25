%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.31 (a)
% Beampattern of Chebychev array with null
% 
% Lillian Xiaolan Xu  9/20/99	
% updated 9/29/00 by K. Bell
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=10; 
lim=-60;
ux_set=-1:0.02:1;
uy_set=-1:0.02:1;

w=[0.773 0.569 0.796 0.029 1;0.569 0.946 0.119 0.618 0.667;0.796 0.119 0.486 0.777 0.286;0.029 0.618 0.777 0.387 0.071;1 0.667 0.286 0.071 0.008];

w=[w(5,:);w(4,:);w(3,:);w(2,:);w(1,:);w];
w=[w(:,5) w(:,4) w(:,3) w(:,2) w(:,1) w];

wd=w(:,1);
for num=2:N
    wd=[wd;w(:,num)];
end

m=(1:N)';

v=exp(i*2*pi*(-m*0.25*0.5-sin(60/180*pi)/2*0.75));
for n=2:N
    v=[v;exp(i*2*pi*(-m*0.25*0.5-n*sin(60/180*pi)/2*0.75))];
end

wo=wd-(wd'*v*inv(v'*v)*v')';

disp('calculating ...')
for x=1:size(ux_set,2)
    ux=ux_set(x);
    for y=1:size(uy_set,2)
        uy=uy_set(y);
        b(x,y)=0;
        for m=1:N
            for n=1:N
                b(x,y)=b(x,y)+wo((n-1)*N+m)*exp(-i*2*pi*(m*ux*0.5+n*uy*0.75));
            end
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
    end
end

figure
mesh(ux_set,uy_set,bdb.');


axis([-1 1 -1 1 lim 0])
%title('Beampatterns of Rectangular Array')
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14);
grid on
