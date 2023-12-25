%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.22
% Beampattern of Circular Array of 20*20, psi_r=0.4*pi;
% Xiaomin Lu 11/2/98	
% Updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all

N = 20;
psi_r = 0.4*pi;
u=-1:1/101:1;
i = 1;
for phi=[0 pi/6 pi/4]

ux = u*cos(phi);
uy = u*sin(phi);


Beam = 0;
for n = 1:N/2
   for m = 1:N/2
      s1 = n-0.5;
      s2 = m-0.5;
      t1 = sqrt(s1^2+s2^2);
       if (t1<=N/2)
          w = psi_r*besselj(1,psi_r*t1)/2/pi/t1;
                   Beam = Beam + 4*w*cos(s1*pi*ux).*cos(s2*pi*uy);
       end
   end
     
end
 
Beam = 20*log10( abs(Beam)/max(abs(Beam)) );
table(i,:) = Beam;
i = i+1;
end

plot(u,table(1,:),'-',u,table(2,:),'--',u,table(3,:),'-.')
h=legend('\phi = 0','\phi = \pi/6','\phi = \pi/4'); 
set(h,'Fontsize',12)
axis([-1 1 -60 0]);
grid
%title('Plot Cut of planar array with circular boundary, N=M=20, \psi_r=0.4*\pi, d=\lambda/2') 
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);




         
      
      
      
 
