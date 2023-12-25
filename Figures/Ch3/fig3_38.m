%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.38
% Comparison of various beam pattern synthesis techniques
% Xin Zhang 3/24/99
% Last updated by K. Bell 9/22/00, 9/30/01
% Functions called: sinc, remez from SP toolbox
%               remezf , remefrf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 11;
M = (N-1)/2;
n = (-M:M).';
u = -1:0.001:1;
vv = exp(j*pi*n*u);
lu = length(u);

%%%%%%%%%%%%%% Define the desired pattern
Bd = zeros(1,lu);
for l = 1:lu
   if ( (acos(u(l)) >= pi/3) & (acos(u(l)) <= 2*pi/3) )
      Bd(l) = 1;
   end
end

%%%%%%%%%%%%%% Parks-McClellan algorithm
F = [0 0.5-0.1 0.5+0.1 1];
A = [1 1 0 0];
W = [1 1];
wr = remez(N-1,F,A,W);
w1 = wr/sum(wr);
B_pm = real(wr*vv);

%%%%%%%%%%%%%% Woodward sampling technique
% Sample the desired pattern
for m = 0:N-1
   ul(m+1) = (m-(N-1)/2)*2/N;        % Equ.(3.118)
   r = abs(ul(m+1));
   if r < 0.5
      Bpsi(m+1) = 1;
   elseif r == 0.5
      Bpsi(m+1) = 0.5;
   else 
      Bpsi(m+1) = 0;
   end
end
Bul = 0;
for m = 0:N-1
   Bul = Bul+Bpsi(m+1)*sinc(N/2*(u-ul(m+1)))./sinc(1/2*(u-ul(m+1)));
end
B_w = Bul;

%%%%%%%%%%%%%% Fourier series method
ad = 0.5*sinc(0.5*n);
w3 = ad/sum(ad);
B_f = real(ad'*vv);

figure
plot(u,B_pm,'-',u,B_w,'--',u,B_f,'-.',ul,Bpsi,'o')
hold on;
plot(u,Bd,':')
hold off;
h=legend('Parks-McClellan','Woodward','Fourier');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)

%axis([-1 1 0 1.201])

figure
plot(u,B_pm,'-',u,B_w,'--',u,B_f,'-.',ul,Bpsi,'o')
hold on;
plot(u,Bd,':')
hold off;
h=legend('Parks-McClellan','Woodward','Fourier');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)
axis([-0.5 0.5 0.9 1.1])