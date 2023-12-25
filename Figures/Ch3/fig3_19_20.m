%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.19
% Figure 3.20
% Array and aperture beam patterns using
% wavenumber sampling technique: symmetric sampling
% Xin Zhang 3/18/99
% Lillian Xu 04/16/2001, K. Bell 7/22/01, 9/30/01
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

u = -1:0.001:1;
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
    if r > 0.5 +eps                % sample the desired pattern
        Bum(m+1) = 0;				
    elseif r > 0.5-eps & r<0.5+eps
        Bum(m+1) = 0.5;
    else 
        Bum(m+1) = 1;
    end
end
Buc = 0;
for m = 0:Ns-1
    Bi = Bum(m+1)*sinc(Ll*(u-up(m+1)));    			% individual pattern
    Buc = Buc+Bi;												% composite pattern
end
%Buc = abs(Buc);

%%%%%%%%%%% Array, sample symmetrically
N = 10;
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
%Bul = abs(Bul);

figure
plot(u,Buc,'-',u,Bul,'--',u,Bd,':',ul,Bpsi,'o')
%title([ 'L = 5.0lambda, N = ' num2str(N) ', d = lambda/2' ])
%axis([-1 1 -0 1.201])
h=legend('Aperture','Array','Desired','Samples');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)

%%%%%%%%%%% Array, sample symmetrically
N = 11;
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
%Bul = abs(Bul);

figure
plot(u,Buc,'-',u,Bul,'--',u,Bd,':',ul,Bpsi,'o')
%title([ 'L = 5.0lambda, N = ' num2str(N) ', d = lambda/2' ])
%axis([-1 1 0 1.201])
h=legend('Aperture','Array','Desired','Samples');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)
