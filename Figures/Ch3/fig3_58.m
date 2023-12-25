%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.58
% Beam pattern:
% N=11, K=8, f1 and f7 plus endpoint references fl and fu
% Lillian Xiaolan Xu 4/13/99
% Updated by K. Bell 7/23/01, 9/30/01
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
f_set=[f1 f1+(f2-f1)/2/K f2-(f2-f1)/2/K f2];
c=335.28;
lam2=c/f2;
N=11;
Na=1; % The number of harmonic nesting structure
d(1)=lam2/2;
n=[-(N-1)/2:(N-1)/2]';
%-------------------------------------------
%-------------------------------------------
for na=1:Na
    for nf=1:4
        f=f_set(na,nf);
        un=c*n/(f*d(na)*N);
        B=sinc(N*un/4)./sinc(un/4);%./sinc(un/4); % B is sinc(Nd/lambda(lower)*un)
        B_con=sinc(N*u/4)./sinc(u/4);
        w=0;
        for k=-(N-1)/2:(N-1)/2
            w=w+B(k+(N-1)/2+1)*exp(i*k*n*2*pi/N)/N;   
        end               % FIB weightings
        w=ones(N,1)/N;   % conventional weightings
        Beam((na-1)*(K+1)+nf,[1:size(u,2)])=w'*exp(i*2*pi*f/c*d(na)*n*u);
    end
end
figure
subplot(2,2,1)
plot(u,20*log10(abs(Beam(1,:))));
grid on
axis([-2 2 -60 0])
title('\it f_l','Fontsize',14)
%xlabel('\it u','Fontsize', 14)
ylabel('Beam pattern (dB)','Fontsize',14)

subplot(2,2,2)
plot(u,20*log10(abs(Beam(2,:))));
axis([-2 2 -60 0])
grid on
title('\it f_0','Fontsize',14)
%xlabel('\it u','Fontsize', 14)
ylabel('Beam pattern (dB)','Fontsize',14)

subplot(2,2,3)
plot(u,20*log10(abs(Beam(3,:))));
grid on
axis([-2 2 -60 0])
title('\it f_7','Fontsize',14)
xlabel('\it u','Fontsize', 14)
ylabel('Beam pattern (dB)','Fontsize',14)

subplot(2,2,4)
plot(u,20*log10(abs(Beam(4,:))));
axis([-2 2 -60 0])
grid on
title('\it f_u','Fontsize',14)
xlabel('\it u','Fontsize', 14)
ylabel('Beam pattern (dB)','Fontsize',14)
