%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 3.29(b)																	
% Taylor & Chebychev  beampattern										
% Xiaomin Lu 11/2/98	
% Updated by K. Bell 9/5/00
% Lillian Xu 04/16/2001, K. Bell 7/29/01,9/30/01
% Functions called: cheby	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

u = -1:1/200:1;

N = 21;
d = 0.5;
scale = N*d;
SLL = 30;
nBar = 6;
A = 1/pi*acosh(10^(SLL/20));

R = 10^(SLL/20);   
x0 = cosh( 1/(N-1)*acosh(R));
x = x0*cos(pi*d*u);
beam = cheby(N-1,x);
beam = abs(beam)/max(abs(beam));
beam = 20*log10(beam);
plot(u,beam,'-')
hold on
 

% Taylor
uroot=[1:1:floor((N-1)/2)]/scale;     % (N-1)/2 for odd and N/2 - 1 for N even
for n = 1:nBar-1
    Vn = nBar*sqrt( (A^2+(n-0.5)^2)/(A^2+(nBar-0.5)^2) );
    uroot(n) = Vn/scale;
end
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

h=legend('Dolph-Chebychev','Taylor');
grid
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
text(-0.03,-78,'(b)','Fontsize',14)









