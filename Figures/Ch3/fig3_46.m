%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.46
% Beam pattern for 4-element MRLA with uniform weighting
% Xin Zhang 3/29/99
% Last updated by K. Bell 9/5/00, 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 4;
Na = 6;
Nan = Na+1;

% 7-element array
n7 = (-(Nan-1)/2:(Nan-1)/2).';
psi = pi*(-1:0.001:1);
v7 = exp(j*n7*psi);
w7 = 1/Nan*[ones(Nan,1)];
B7 = w7'*v7;
B7 = 10*log10(abs(B7).^2);

% MRLA
nm4 = [-3 -2 1 3]';
v4 = exp(j*nm4*psi);
w4 = 1/N*[ones(N,1)];
B4 = w4'*v4;
B4 = 10*log10(abs(B4).^2);

% Output the results
figure
plot(psi/pi,B4,'-')
hold on;
plot(psi/pi,B7,'-.')
hold off;
axis([-1 1 -40 0])
grid on;
h=legend('MRLA(132), {\it N}=4','Conventional, {\it N_a}=7');
set(h,'Fontsize',14)
xlabel('\it u','Fontsize',16)
ylabel('Beam pattern (dB)','Fontsize',16)

% Find HPBW for the non-uniform 4-element array
u = [0:0.0001:1];
vv = exp(j*nm4*pi*u);
B = w4'*vv;
g = find(B<1/sqrt(2));
HPBW_u_m4 = 2*u(min(g));
HPBW_psi_m4 = pi*HPBW_u_m4;
% Find HPBW for the uniform 4-element array
n4 = (-(N-1)/2:(N-1)/2).';
vv = exp(j*n4*pi*u);
B = w4'*vv;
g = find(B<1/sqrt(2));
HPBW_u_4 = 2*u(min(g));
HPBW_psi_4 = pi*HPBW_u_4;
% Find HPBW for the uniform 7-element array
vv = exp(j*n7*pi*u);
B = w7'*vv;
g = find(B<1/sqrt(2));
HPBW_u_7 = 2*u(min(g));
HPBW_psi_7 = pi*HPBW_u_7;

% Find BWmm for the non-uniform 4-element array
u = [0.2:0.0001:0.3];
vv = exp(j*nm4*pi*u);
B = w4'*vv;
[y,I]=min(abs(B));
umin = u(I);
BWmm_u_m4 = 2*umin;
psimin = pi*umin;
BWmm_psi_m4 = 2*psimin;
Bmin_m4 = y;
% Find BWnn for the uniform 4-element array
u = [0:0.0001:1];
vv = exp(j*n4*pi*u);
B = w4'*vv;
g = find(B<0);
null = u(min(g));
BWNN_psi_4 = 2*null*pi;
% Find BWnn for the uniform 7-element array
vv = exp(j*n7*pi*u);
B = w7'*vv;
g = find(B<0);
null = u(min(g));
BWNN_psi_7 = 2*null*pi;
