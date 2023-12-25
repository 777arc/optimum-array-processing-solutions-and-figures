%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 3.29(a)																	
% Continuous & Sampled Taylor pattern									
% Xiaomin Lu 11/2/98	
% Updated by K. Bell 9/5/00
% Lillian Xu 04/16/2001, K. Bell 7/31/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

N = 21;

u = -1:1/200:1;
d = 0.5;
scale = N*d;
SLL = 30;
nBar = 6;

v = scale*u;
A = 1/pi*acosh(10^(SLL/20));
uroot=[1:1:floor((N-1)/2)]/scale;     % (N-1)/2 for odd and N/2 - 1 for N even

Bv = sinc(v);
for n = 1:nBar-1
    Vn = nBar*sqrt( (A^2+(n-0.5)^2)/(A^2+(nBar-0.5)^2) );
    uroot(n) = Vn/scale;
    Bv = Bv.*(1-v.^2/Vn^2)./(1-v.^2/n^2);
end
Bv = abs(Bv)/max(abs(Bv));  									% Continuous pattern
plot(u,20*log10(Bv),'-')


hold on

if rem(N,2) == 0 % even
    roots = [0 uroot -uroot 1];
else
    roots = [0 uroot -uroot];
end

n = (-(N-1)/2:(N-1)/2)';
Vr = exp(j*n*pi*roots);

w = inv(Vr')*[1;zeros(N-1,1)];

beam = w'*exp(j*n*pi*u);
beam = abs(beam)/max(abs(beam));
beam = 20*log10(beam);
plot(u,beam,'--')													%zero matching


grid
%title('Continuous and sampled Taylor pattern, N=21, nBar=6, d=\lambda/2, SLL=-30dB');
h=legend('Aperture','Array');%,'Zero matching');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
text(-0.03,-78,'(a)','Fontsize',14)












