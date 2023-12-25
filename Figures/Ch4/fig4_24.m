%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4.24
% Beampattern of Circular Array of 20*20, psi_r=0.4*pi, Kaiser Window
% Xiaomin Lu 11/2/98	
% updated by K. Bell 9/29/00
% Updated by Lillian Xiaolan Xu 02/12/2001, K. Bell 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 20;
beta = 5;
psi_r = 0.4*pi;
u=-1:1/101:1;
i = 1;
for phi=0
ux = u*cos(phi);
uy = u*sin(phi);


Beam1 = 0;
Beam2 = 0;
for n = 1:N/2
   for m = 1:N/2
      s1 = n-0.5;
      s2 = m-0.5;
      t1 = sqrt(s1^2+s2^2);
       if (t1<=N/2)
          w = psi_r*besselj(1,psi_r*t1)/2/pi/t1;
          Beam1 = Beam1 + 4*w*cos(s1*pi*ux).*cos(s2*pi*uy);
          w = w*besseli(0,beta*sqrt(1-(t1/14)^2))/besseli(0,beta);
          Beam2 = Beam2 + 4*w*cos(s1*pi*ux).*cos(s2*pi*uy);
       end
   end
end
 
Beam1 = 20*log10( abs(Beam1)/max(abs(Beam1)) );
Beam2 = 20*log10( abs(Beam2)/max(abs(Beam2)) );
end

%plot(u,Beam1,'-',u,Beam2,'--')
%legend('Uniform','Kaiser Win'); 
plot(u,Beam2);
axis([-1 1 -60 0]);
grid
%title('Beam of planar array, Kaiser Window, N=M=20, \psi_r=0.4*\pi, d=\lambda/2') 
xlabel('\it u_r','Fontsize',14);
ylabel('Beam pattern (dB)','Fontsize',14);





         
      
      
      
 
