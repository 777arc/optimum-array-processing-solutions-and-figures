%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.13
% Xiaomin Lu
% Updated by K. Bell 10/13/00, 7/23/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distribution of eigenvalues

clear all
close all
N = 10;

du = [0.0433 0.5];

for n=1:2
   phaseRange = pi*[0 0.5 0.75 1];
   rhoRange = 0:1/500:0.999;
   
   %B12 = diric(du(n)*pi, N);
   B12 = sin(N*pi*du(n)/2)./(N*sin(pi*du(n)/2));

   S1 = 1;
   S2 = S1;
   
   h1 = 1;
   for phase = phaseRange
      h2 = 1;
      for rho = rhoRange
         S12 = sqrt(S1*S2)*rho*exp(j*phase);
         
         b = N*( S1 + S2 + 2*real(S12*B12') );
         c = N^2*S1*S2*(1-rho^2)*(1-B12^2);
         
         g(1,h1,h2) = 0.5*(b + sqrt(b^2-4*c));
         g(2,h1,h2) = 0.5*(b - sqrt(b^2-4*c));
         
         h2 = h2 + 1;
      end
      h1 = h1 + 1;
   end
   
   g = 10*log10(g/N);
   x = rhoRange;
   
   g1 = squeeze(g(1,:,:));
   g2 = squeeze(g(2,:,:));
   figure
   plot(x,g1(1,:),'-',x,g1(2,:),'--',x,g1(3,:),'-.',x,g1(4,:),':');
   hold on
   plot(x,g2(1,:),'-',x,g2(2,:),'--',x,g2(3,:),'-.',x,g2(4,:),':');
   grid
   xlabel('\rho','Fontsize',14);
   ylabel('Normalized eigenvalues (dB)','Fontsize',14);
   %title(['N=10,S1=S2=0dB, \Deltau=' num2str(du(n))])
   h=legend('\phi=0','\phi=\pi/2','\phi=3\pi/4','\phi=\pi',3);
   set(h,'Fontsize',12)
end



