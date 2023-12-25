%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.57
% Nested array
% Figure 3.59
% Desired beamformer
% Figure 3.61
% Beam patterns for broadband array:
% N=11, 500<=f<=8000, constant mainlobe
% 
% Lillian Xu 04/17/2001
% Updated by K. Bell 7/23/01, 9/30/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Harmonic nesting          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

u1=-2:0.01:2;
K=8;

f0=500;
f1=4000;
f2=8000;
f_set(1,:)=[f2:-(f2-f1)/K:f1];
f_set(2,:)=f_set(1,:)/2;
f_set(3,:)=f_set(2,:)/2;
f_set(4,:)=f_set(3,:)/2;
fset=[f_set(1,:) f_set(2,:) f_set(3,:) f_set(4,:)];

c=335.28;

lam2=c/f2;

N=11;
Na=4; % The number of harmonic nesting structure
d(1)=lam2/2;
d(2)=2*d(1);
d(3)=2*d(2);
d(4)=2*d(3);

n=[-(N-1)/2:(N-1)/2]';
figure
for na=1:Na
   plot(n*d(na)/d(1),na*ones(N,1),'*')
   hold on
   plot(n*d(na)/d(1),zeros(N,1),'o')
end
axis([-42 42 -1 5])
%title('harmonic nesting subarrays')

   B_con=sinc(N*u1/4)./sinc(u1/4);
u=-1:0.01:1;
%-------------------------------------------
%-------------------------------------------
for na=1:Na
for nf=1:K+1
   f=f_set(na,nf);
   un=c*n/(f*d(na)*N);
   B=sinc(N*un/4)./sinc(un/4);%./sinc(un/4); % B is sinc(Nd/lambda(lower)*un)
   %  B_con=sinc(N*u/4);  % aperture

   
%  figure
%  plot(u,20*log10(B_con),'--',un,20*log10(B),'o')
%  grid on
   w=0;
   for k=-(N-1)/2:(N-1)/2
      w=w+B(k+(N-1)/2+1)*exp(i*k*n*2*pi/N)/N;   
   end
 %  figure
 %  plot([1:11],w)

   Beam((na-1)*(K+1)+nf,[1:size(u,2)])=w'*exp(i*2*pi*f/c*d(na)*n*u);
%B(nf,:)=B(nf,:)/B(nf,101);   
end
end

  figure
  plot(u1,20*log10(abs(B_con)))
  grid on
xlabel('\it u','Fontsize', 16)
ylabel('Beam pattern (dB)','Fontsize',16)
%title('desired beam pattern')
figure
mesh(u,fset,max(20*log10(abs(Beam)),-60));
axis([-1 1 f0 f2 -60 0])
ylabel('\it f','Fontsize',14)
xlabel('\it u','Fontsize',14)
zlabel('Beam pattern (dB)','Fontsize',14)
%title(['N=',num2str(N),...
%      ', K=',num2str(K)]) %,...
  %    ', c=',num2str(c),'m/s'])