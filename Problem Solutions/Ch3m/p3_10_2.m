%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% problem 3.10.2 
% K. Bell 11/3/00
% Updated 10/10/03 by K. Bell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
N=32;
d=0.5;
D = d*[-(N-1)/2:1:(N-1)/2].';
u = [-1:0.001:1];
v = exp(j*2*pi*D*u);
nu = length(u);
u2 = [-3:1:3]*2/N;
n2 = length(u2);
v2 = exp(j*2*pi*D*u2);


% Taylor weights
nbar = 6;
SL = [-20];

w = zeros(N,1);

R = 10^(-SL/20);
Asq = (acosh(R)/pi)^2;
if rem(N,2)==0  % even
   psi = (2*pi/N)*[1:1:(N/2)-1].';
   psinbar = (2*pi/N)*nbar*sqrt((Asq+[0.5:1:nbar-1.5].^2).'/(Asq+(nbar-0.5).^2));
   psin = [psinbar;psi(nbar:(N/2)-1)];
   psit = [psin;-psin;-pi];
else % odd
   psi = (2*pi/N)*[1:1:(N-1)/2].';
   psinbar = (2*pi/N)*nbar*sqrt((Asq+[0.5:1:nbar-1.5].^2).'/(Asq+(nbar-0.5).^2));
   psin = [psinbar;psi(nbar:(N-1)/2)];
   psit = [psin;-psin];
end
z = exp(j*psit);
wt = poly(z).';

w = wt/sum(wt);

% Seven steered taylor beams
Bno = v2.*(w*ones(1,7));
% check orthonormality
Bno'*Bno

% Hw
C = Bno'*Bno;
[vv,lam] = eig(C);
Hw = vv*inv(sqrt(real(lam)))*vv';

% alternate Hw
%Hw = inv(sqrtm(C));

Bbs = Bno*Hw;
% check orthonormality
Bbs'*Bbs


Bp = Bbs'*v;

% normalize for plotting
c = max(max(abs(Bp)));
Bp = Bp/c;

subplot(2,1,1);
for n=1:7
h=plot(u, 10*log10(abs(Bp(n,:)).^2));
hold on
end
axis([-1 1 -40 0])
ylabel('Beampattern (dB)')
xlabel('u')
title('Problem 3.10.2, Seven Orthogonalized Taylor Beams')
grid on

subplot(2,1,2);
h1=plot(u, 10*log10(abs(Bp(4,:)).^2));
hold on

Bt = Bno(:,4)'*v;
% normalize for plotting
c = max(max(abs(Bt)));
Bt = Bt/c;
h2=plot(u, 10*log10(abs(Bt.^2)),'--');


axis([-1 1 -40 0])
ylabel('Beampattern (dB)')
xlabel('u')
grid on
h = [h2;h1];
legend(h,'Taylor Beam','Orthogonalized Taylor Beam')
