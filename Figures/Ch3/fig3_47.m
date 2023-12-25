%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3.47
% Beam pattern for 5-element MRLA
% Xin Zhang 3/29/99
% Last updated by K. Bell 9/5/00, 7/23/01, 9/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

N = 5;
Na = 9;
Nan = Na+1;

% 10-element array
n10 = (-(Nan-1)/2:(Nan-1)/2).';
psi = pi*(-1:0.001:1);
v10 = exp(j*n10*psi);
w10 = 1/Nan*[ones(Nan,1)];
B10 = w10'*v10;
B10 = 10*log10(abs(B10).^2);

% MRLA(1332)
nm51 = [-4.5 -3.5 -0.5 2.5 4.5]';
v51 = exp(j*nm51*psi);
w5 = 1/N*[ones(N,1)];
B51 = w5'*v51;
B51 = 10*log10(abs(B51).^2);

% MRLA(3411)
nm52 = [-4.5 -1.5 2.5 3.5 4.5]';
v52 = exp(j*nm52*psi);
B52 = w5'*v52;
B52 = 10*log10(abs(B52).^2);

% Output the results
figure
plot(psi/pi,B51,'-')
hold on;
plot(psi/pi,B10,'-.')
hold off;
axis([-1 1 -30 0])
grid on;
h=legend('MRLA(1332), {\it N}=5','Conventional, {\it N_a}=10');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)

figure
plot(psi/pi,B52,'-')
hold on;
plot(psi/pi,B10,'-.')
hold off;
axis([-1 1 -30 0])
grid on;
h=legend('MRLA(3411), {\it N}=5','Conventional, {\it N_a}=10');
set(h,'Fontsize',12)
xlabel('\it u','Fontsize',14)
ylabel('Beam pattern (dB)','Fontsize',14)
