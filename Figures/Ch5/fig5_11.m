%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5.11
% Eigenbeams for two sources
% K. Bell 7/24/01, 10/23/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%************************
% Array
%************************
N = 10;                                 % Elements in array
d = 0.5;                                % sensor spacing half wavelength wrt wc
D = [-(N-1)/2:1:(N-1)/2].';
BWNN = 2/(N*d);
u=[-1:0.001:1];
nu=length(u);
vv = exp(j*pi*D*u);

nf = N;
%************************
% Source
%************************
dus = [2.3 0.68 0.1]*BWNN;

ns = length(dus);
val1 = zeros(1,ns);
val2 = zeros(1,ns);
for n=1:ns
    us = [-dus(n)/2 dus(n)/2];
    AS = exp(j*2*pi*d*D*us)/sqrt(N);
    
    v1 = AS(:,1)+AS(:,2);      % sum
    v1 = v1/norm(v1);          % normalized
    v2 = AS(:,1)-AS(:,2);      % difference
    v2 = v2/norm(v2);          % normalized

    B1 = real(v1'*vv);
    B2 = real(v2'*vv);
    figure
    plot(u,B1/10,'-');
    hold on
    plot(u,B2/10,'--');
    
    plot(us(1)*[1 1],[-2 2],'-g')
    plot(us(2)*[1 1],[-2 2],'-g')
    xlabel('\itu','Fontsize',14)
    ylabel('Eigenbeams','Fontsize',14)
    grid on
    hold off
    axis([-1 1 -0.4 0.4])
    set(gca,'XTick',[-1:0.2:1])
    
end
