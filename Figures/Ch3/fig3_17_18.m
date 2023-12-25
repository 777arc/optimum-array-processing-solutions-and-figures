%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.17
% Pattern sampling
% Figure 3.18
% Desired and synthesized patterns
% Woodward synthesis procedure
% Xin Zhang 3/18/99
% Updated 9/14/00, K. Bell 7/22/01, 9/30/01
% Functions called: sinc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

Ll = 5.0;
u = -1:0.001:1;
lu = length(u);

% Define the desired pattern
Bd = zeros(1,lu);
for l = 1:lu
    if ( (acos(u(l)) >= pi/3) & (acos(u(l)) <= 2*pi/3) )
        Bd(l) = 1;
    end
end

% Sample the desired pattern
delta_us = 1/Ll;     			% sample interval
Ns = floor(2*Ll);        		% number of samples
for m = 0:Ns-1
    up(m+1) = m*delta_us-1;		% values of u at sample points
end
for m = 0:Ns-1						% Find the values of pattern at sample points
    r = abs(up(m+1));
    if r > 0.5
        Bum(m+1) = 0;
    elseif r == 0.5
        Bum(m+1) = 0.5;
    else 
        Bum(m+1) = 1;
    end
end

figure
plot(u,Bd,'-',up,Bum,'o')
xlabel('\it u','Fontsize',14)
ylabel('\it B(u)','Fontsize',14)
axis([-1.2 1.2 0 1.2])

figure
h1 = plot(u,Bd,':');
hold on;

% Woodward synthesis, 
Bu = 0;
for m = 0:Ns-1
    Bi = Bum(m+1)*sinc(Ll*(u-up(m+1)));    		% individual pattern
    Bu = Bu+Bi;												% composite pattern
    h2 = plot(u,Bi,'-');
end
h3 = plot(u,Bu,'--');
hold off;
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern','Fontsize',14)
h=legend([h2 h3 h1],'Individual','Composite','Desired');
set(h,'Fontsize',12)