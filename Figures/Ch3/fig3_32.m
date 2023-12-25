%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.32
% Desired beam pattern and synthesized patterns 
% using the Fourier series method and Woodward 
% sampling technique
% (a) N = 10; (b) N = 20
% Xin Zhang 3/23/99
% K. Bell 9/14/00, K. Bell 7/22/01, 9/30/01
% Function called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

psi0 = 0.5*pi;
psi = -1*pi:0.001*pi:pi;
u = -1:0.001:1;
lu = length(u);

%%%%%%%%%%% Define the desired pattern
Bd = zeros(1,lu);
for l = 1:lu
   if ( (acos(u(l)) >= pi/3) & (acos(u(l)) <= 2*pi/3) )
      Bd(l) = 1;
   end
end

N = 10;
%%%%%%%%%%%%%%% Fourier series method
for m = (-N/2+1):N/2
   a(m+N/2) = psi0/pi*sinc((m-1/2)*psi0/pi);
end

s = 0;
for m = 1:N/2
   s = s+a(m+N/2)*cos((m-1/2)*psi);
end
Bt = 2*s;                      % truncated beam pattern
Bta = abs(Bt);

%%%%%%%%%%%%%% Woodward sampling
% Sample the desired pattern
for m = 0:N-1
   ul(m+1) = (m-(N-1)/2)*2/N;        
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
Bwa = abs(Bul);

figure
plot(psi/pi,Bta,'-',psi/pi,Bwa,'--',psi/pi,Bd,':')
h=legend('Fourier series','Woodward sampling','Desired');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('|{\it B(u)}|','Fontsize',14)

N = 20;
%%%%%%%%%%%%%%% Fourier series method
for m = (-N/2+1):N/2
   a(m+N/2) = psi0/pi*sinc((m-1/2)*psi0/pi);
end
%a/max(a)

s = 0;
for m = 1:N/2
   s = s+a(m+N/2)*cos((m-1/2)*psi);
end
Bt = 2*s;                      % truncated beam pattern
Bta = abs(Bt);

%%%%%%%%%%%%%% Woodward sampling
% Sample the desired pattern
for m = 0:N-1
   ul(m+1) = (m-(N-1)/2)*2/N;        
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
Bwa = abs(Bul);

figure
plot(psi/pi,Bta,'-',psi/pi,Bwa,'--',psi/pi,Bd,':')
h=legend('Fourier series','Woodward sampling','Desired');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('|{\it B(u)}|','Fontsize',14)
