% Script problem6613.m: Solves Problem 6.6.13 parts (a) and (c) in Van Trees, Volume IV.
% 11/16/03, J.A. Tague

% Initialize parameters for a ten-element standard line array

N = 10;
v_s = ones(10, 1);
ind = 0 : N - 1;
u = -1 : 1/512 : 1 - 1/512;

% Find the array's broadband noise cross-spectral density matrix

B_fI = 0.4; % The normalized temporal interference bandwidth
u_I1 = 0.3; % DOA from first interference source
u_I2 = 0.5; % DOA from second interference source

row1 = exp(j * pi * ind * u_I1) .* sinc(0.5 * ind * B_fI * u_I1);
S_1 = toeplitz(row1);

row1 = exp(j * pi * ind * u_I2) .* sinc(0.5 * ind * B_fI * u_I2);
S_2 = toeplitz(row1);

S_ntrue = 100 * S_1 + 100 * S_2 + eye(10);
rho_n = N * S_ntrue / trace(S_ntrue);

% Find the optimum MVDR beamformer and calculate array gain, assuming that the noise is 
% narrowband

n = N - 1 : -1 : 0;
n = n';

v1 = exp(j * pi * u_I1 * (n - (N - 1) / 2));
v2 = exp(j * pi * u_I2 * (n - (N - 1) / 2));
S_ndesign = 100 * v1 * v1' + 100 * v2 * v2' + eye(10);

w = inv(S_ndesign) * v_s / (v_s' * inv(S_ndesign) * v_s);

AG = abs(w' * v_s)^2 / (w' * rho_n * w);
AG = 10 * log10(real(AG))

% Plot the beam pattern of the MVDR beamformer, designed assuming that the
% interference sources are narrowband

W = fft(conj(w), 1024);
W = fftshift(W);
B = 20 * log10(abs(W));

for k1 = 1 : size(B)
    if (B(k1) < -90) B(k1) = -90;
    end;
end;

figure;
plot(u, B);
xlabel('Direction Cosine u')
ylabel('Beam Pattern (dB)')
title('Beam Pattern of MVDR Beamformer Using Narrowband Noise Covariance Matrix: u_I = 0.3 and 0.5');
grid;

% Repeat using the actual noise spectral density matrix to determine w

w = inv(S_ntrue) * v_s / (v_s' * inv(S_ntrue) * v_s);

% Calculate the array gain this beamformer provides.  It should be more
% than that of the narrowband assumption if B_fI > 0.

AG = abs(w' * v_s)^2 / (w' * rho_n * w);
AG = 10 * log10(real(AG))

% Plot the beam pattern of the MVDR beamformer that uses the actual noise
% spectral density matrix to determine w

W = fft(conj(w), 1024);
W = fftshift(W);
B = 20 * log10(abs(W));

for k1 = 1 : size(B)
    if (B(k1) < -90) B(k1) = -90;
    end;
end;

figure;
plot(u, B);
xlabel('Direction Cosine u')
ylabel('Beam Pattern (dB)')
title('Beam Pattern of MVDR Beamformer Using Actual Noise Covariance Matrix: u_I = 0.3 and 0.5; B_{fI} = 0.4');
grid;