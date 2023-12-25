%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.60
% Samples and array weights; N=11, K=8, conventional desired pattern 
% Lillian Xiaolan Xu 4/13/99
% last updated by K. Bell 7/23/01, 9/30/01, 10/17/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Harmonic nesting          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

u=-2:0.01:2;
K=8;

f1=4000;
f2=8000;

f_set=[f1+(f2-f1)/2/K f2-(f2-f1)/2/K f2];

c=335.28;

lam2=c/f2;

N=11;
Na=1; % The number of harmonic nesting structure
d(1)=lam2/2;

n=[-(N-1)/2:(N-1)/2]';

%-------------------------------------------
%-------------------------------------------
figure
for na=1:Na
   for nf=1:2
      f=f_set(na,nf);
      un=c*n/(f*d(na)*N);
      B=sinc(N*un/4)./sinc(un/4);%./sinc(un/4); % B is sinc(Nd/lambda(lower)*un)
      B_con=sinc(N*u/4)./sinc(u/4);
      
      subplot(2,2,1+(nf-1)*2)
      plot(u,20*log10(abs(B_con)),'--',un,max(20*log10(abs(B)),-60),'o')
      grid on
      ylabel('Beam pattern (dB)','Fontsize',14)
      axis([-2 2 -60 0])
      
      w=0;
      for k=-(N-1)/2:(N-1)/2
         w=w+B(k+(N-1)/2+1)*exp(i*k*n*2*pi/N)/N;   
      end               % FIB weightings
      
      subplot(2,2,2+(nf-1)*2) 
      plot([1:11],real(w))
      axis([1 11 -0.05 0.2])
      ylabel('Weights','Fontsize',14)
      
      Beam((na-1)*(K+1)+nf,[1:size(u,2)])=w'*exp(i*2*pi*f/c*d(na)*n*u);
   end
end

subplot(2,2,1)
title('Samples','Fontsize',14)
subplot(2,2,3)
xlabel('\it u','Fontsize',14)
subplot(2,2,2)
title('Weights','Fontsize',14)
text(11.3,0.075,'f_0','Fontsize',14)
subplot(2,2,4)
xlabel('Sensors','Fontsize',14)
text(11.3,0.075,'f_7','Fontsize',14)

