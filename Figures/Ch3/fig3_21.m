%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.21
% Array and aperture beam patterns using
% wavenumber sampling technique: 
% d = lambda/4, symmetric samples
% Xin Zhang, 3/21/99
% Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% d = lambda/4 => differences:
% 1) u = -2:0.001:2
% 2) ul(m+1) = (m-(N-1)/2)*4/N
% 3) if use direct imposing, then 
%       Bul = Bul+Bpsi(m+1)*sinc(N/4*(u-ul(m+1)))./sinc(1/4*(u-ul(m+1)))
%    if use DFT, then
%       Vo = exp(j*amf*pi*u/2); while psik does not change

u = -2:0.001:2;
lu = length(u);

%%%%%%%%%%% Define the desired pattern
Bd = zeros(1,lu);
for l = 1:lu
    if ( (acos(u(l)) >= pi/3) & (acos(u(l)) <= 2*pi/3) )
        Bd(l) = 1;
    end
end

%%%%%%%%%%% Aperture
Ll = 5.0;
ofs = 0.5/Ll;
eps = 0.001;
% Sample the desired pattern
delta_us = 1/Ll;     			% sample interval
Ns = floor(2*Ll);        		% number of samples
for m = 0:Ns-1
    up(m+1) = m*delta_us-1+ofs;		% values of u at sample points
    r = abs(up(m+1));
    %up(m+1) = m*delta_us-1;		% values of u at sample points
    %r = abs(up(m+1));
    if r > 0.5 +eps                % sample the desired pattern
        Bum(m+1) = 0;				
    elseif r > 0.5-eps & r<0.5+eps
        Bum(m+1) = 0.5;
    else 
        Bum(m+1) = 1;
    end
    %   if r > 0.5                 % sample the desired pattern
    %      Bum(m+1) = 0;				
    %   elseif r == 0.5
    %      Bum(m+1) = 0.5;
    %   else 
    %      Bum(m+1) = 1;
    %   end
end
Buc = 0;
for m = 0:Ns-1
    Bi = Bum(m+1)*sinc(Ll*(u-up(m+1)));    			% individual pattern
    Buc = Buc+Bi;												% composite pattern
end

%%%%%%%%%%% Array, sample symmetrically
N = 20;
% Sample the desired pattern
for m = 0:N-1
    ul(m+1) = (m-(N-1)/2)*4/N;
    r = abs(ul(m+1));
    % sample the desired beam pattern
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
    Bul = Bul+Bpsi(m+1)*sinc(N/4*(u-ul(m+1)))./sinc(1/4*(u-ul(m+1)));
end
% or
for k = 0:N-1
    psik = (k-(N-1)/2)*2*pi/N;                   % psik's do not change
    B(k+1) = exp(-j*(N-1)/2*psik)*Bpsi(k+1)';    %   because they are still from -pi to pi
end

for n = 0:N-1
    b(n+1) = 0;
    for k = 0:N-1
        b(n+1) = b(n+1)+B(k+1)*exp(j*k*n*2*pi/N);
    end
end
b = 1/N*b;                       % b = ifft(B);

for n = 0:N-1
    w(n+1) = b(n+1)*exp(-j*pi/N*(N-1)*n);
end 
w=w';

amf = [-(N-1)/2:(N-1)/2]';
Vo = exp(j*amf*pi*u/2);
Beaml = w'*Vo;

figure
plot(u,Buc,'-',u,Bul,'--',u,Bd,':',ul,Bpsi,'o')
%title('L = 5.0lambda, N = 20, d = lambda/4')
h=legend('Aperture','Array','Desired','Samples');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)
