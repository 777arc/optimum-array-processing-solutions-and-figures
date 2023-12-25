% Script problem6610b.m: Solves Problem 6.6.10 part (b) in Van Trees, Volume IV
% 11/16/03, J.A. Tague

% Initialize parameters for a ten element standard line array

N = 10;
n = N - 1 : -1 : 0;
n = n';

% Assume the signal arrives from 30 degrees and calculate its array
% manifold vector.  That vector will be used to build the optimum
% beamformer.

theta = 30 * pi / 180;
v_m = exp(j * pi * cos(theta) * (n - (N - 1) / 2));

SNR_db = -30 : 1 : 30;
AG_1 = zeros(1, length(SNR_db));
AG_2 = zeros(1, length(SNR_db));

% Calculate array gain assuming f = 1.2 * f_c (mismatched frequency)

delta = 0.2;
v_a = exp(j * pi * cos(theta) * (1 + delta) * (n - (N - 1) / 2));

for k = 1 : length(SNR_db)
    SNR = 10^(SNR_db(k) / 10);
    S_x = SNR * v_a * v_a' + eye(N);
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG_1(k) = abs(w' * v_a)^2 / (w' * w);
    AG_1(k) = 10 * log10(real(AG_1(k)));
end

% Calculate array gain assuming f = 1.4 * f_c (mismatched frequency)

delta = 0.4;
v_a = exp(j * pi * cos(theta) * (1 + delta) * (n - (N - 1) / 2));

for k = 1 : length(SNR_db)
    SNR = 10^(SNR_db(k) / 10);
    S_x = SNR * v_a * v_a' + eye(N);
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG_2(k) = abs(w' * v_a)^2 / (w' * w);
    AG_2(k) = 10 * log10(real(AG_2(k)));
end

% Plot the results

figure;
plot(SNR_db, AG_1, SNR_db, AG_2);
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Array Gain (dB)');
title('Array Gain: MPDR Beamformer with Frequency Mismatches');
legend('f = 1.2 f_c', 'f = 1.4 f_c');

