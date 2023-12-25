%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.31 (b)
% Beampattern of Chebychev array with null
% Lillian Xiaolan Xu 9/20/99	
% Updated by K. Bell 9/29/00, 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N=10; 
lim=-60;
ux_set=-1:0.02:1;
uy_set=-1:0.02:1;
theta_set=pi*[-0.5:0.002:0.5];
phi_set=[55 60 65];

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
fig=1;
figure
for num=1:size(phi_set,2)
   phi=phi_set(num);
   for x=1:size(theta_set,2)
      theta=theta_set(x);
      ux=sin(theta)*cos(phi/180*pi);
      uy=sin(theta)*sin(phi/180*pi);      
      b(x)=0;
      for m=1:N
         for n=1:N
            b(x)=b(x)+wo((n-1)*N+m)*exp(-i*2*pi*(m*ux*0.5+n*uy*0.75));
         end
      end
   end
   b=b/max(max(b));
   bdb=20*log10(abs(b));
   
   for x=1:size(ux_set,2)
      if bdb(x)<lim
         bdb(x)=lim;
      end
   end
   
   subplot(3,1,fig);
   fig=fig+1;
   plot(theta_set/pi*180,bdb);
   axis([-90 90 lim 0])
   title(['\phi=-',num2str(phi),'^o'],'Fontsize',14)
   ylabel('Beam pattern (dB)','Fontsize',14)
   if fig==4
      xlabel(['\theta^o'],'Fontsize',14);
   end
   grid on
end;
