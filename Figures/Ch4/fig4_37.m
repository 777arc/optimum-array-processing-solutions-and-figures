%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.37
% Beampatterns for a uniform circular array of 20 elements
% Lillian Xiaolan Xu 9/20/99	
% Last updated by K. Bell 9/30/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01, 11/15/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
N=20; 
lim=-20;
ux_set=-1:0.02:1;
uy_set=-1:0.02:1;
nx = size(ux_set,2);
ny = size(uy_set,2);
b = zeros(nx,ny);

disp('calculating...')
for x=1:nx
   ux=ux_set(x);
   for y=1:ny
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
h=mesh(uy_set,ux_set,bdb);
axis([-1 1 -1 1 lim 0])
%title('Beampatterns of Circular Array')
xlabel('\it u_x','Fontsize',14);
ylabel('\it u_y','Fontsize',14);
zlabel('Beam pattern (dB)','Fontsize',14);
grid on
colormap([0 0 0])