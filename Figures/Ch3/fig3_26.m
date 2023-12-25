%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.26
% Patterns of Dolph-Chebychev arrays with eight elements
%      (-20, -30, -40 dB sidelobes)
% Lillian Xiaolan Xu 3/23/99
% updated by K. Bell 7/22/01, 9/30/01, 10/17/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%   for N=8
N=8;
sidelobe=[20 30 40];   	% *dB below the main lobe maximum
R=10.^(sidelobe/20);
x0_set=cosh(1/(N-1)*acosh(R));

u = -1:1/200:1;

for num=1:size(x0_set,2)
    x0=x0_set(num);
    a3=x0^7;
    a2=(a3-x0^5)*112/16;
    a1=((x0^3-a3)*56+20*a2)/4;
    a0=-7*x0+3*a1-5*a2+7*a3;
    
    wdq=[a3 a2 a1 a0 a0 a1 a2 a3];        %w'
    sumw=sum(wdq);
    wdq=wdq./sumw;
    wdq=wdq';
    b=0;     
    for m=1:N
        b=b+wdq(m)*exp(i*(-(N+1)/2+m)*pi*u);
    end;
    bdb1(num,:)=20*log10(abs(b));
    
    W(:,num)=wdq;
    
    %x=x0*cos(pi*u/2);
    %T7=64*x.^7-112*x.^5+56*x.^3-7*x;
end;

figure
plot(u,bdb1(1,:),'-',u,bdb1(2,:),'--',u,bdb1(3,:),'-.');
axis([-1 1 -60 10])
h=legend('Sidelobes=-20 dB','Sidelobes=-30 dB','Sidelobes=-40 dB');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
grid on