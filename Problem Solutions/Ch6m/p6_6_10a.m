% Script problem6610a.m: Solves Problem 6.6.10 part (a) in Van Trees, Volume IV
% 11/16/03, J.A. Tague

% Initialize parameters for a ten element standard line array

N = 10;
n = N - 1 : -1 : 0;
n = n';

delta = -1 : 0.01 : 1;

% Assume the signal arrives from 30 degrees and calculate its array
% manifold vector.  That vector will be used to build the optimum
% beamformer.

theta = 30 * pi / 180;
v_m = exp(j * pi * cos(theta) * (n - (N - 1) / 2));

AG_1 = zeros(1, length(delta));
AG_2 = zeros(1, length(delta));
AG_3 = zeros(1, length(delta));
AG_4 = zeros(1, length(delta));
AG_5 = zeros(1, length(delta));

% Calculate array gain for SNR = -10 dB

for k = 1 : length(delta)
    SNR = 10^(-10 / 10);
    v_a = exp(j * pi * cos(theta) * (1 + delta(k)) * (n - (N - 1) / 2));
    S_x = SNR * v_a * v_a' + eye(N);
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG_1(k) = abs(w' * v_a)^2 / (w' * w);
    AG_1(k) = 10 * log10(real(AG_1(k)));
end

% Calculate array gain for an SNR = 0 dB

for k = 1 : length(delta)
    SNR = 10^(0 / 10);
    v_a = exp(j * pi * cos(theta) * (1 + delta(k)) * (n - (N - 1) / 2));
    S_x = SNR * v_a * v_a' + eye(N);
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG_2(k) = abs(w' * v_a)^2 / (w' * w);
    AG_2(k) = 10 * log10(real(AG_2(k)));
end

% Calculate array gain for SNR = 10 dB

for k = 1 : length(delta)
    SNR = 10^(10 / 10);
    v_a = exp(j * pi * cos(theta) * (1 + delta(k)) * (n - (N - 1) / 2));
    S_x = SNR * v_a * v_a' + eye(N);
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG_3(k) = abs(w' * v_a)^2 / (w' * w);
    AG_3(k) = 10 * log10(real(AG_3(k)));
end

% Calculate array gain for an SNR = 20 dB

for k = 1 : length(delta)
    SNR = 10^(20 / 10);
    v_a = exp(j * pi * cos(theta) * (1 + delta(k)) * (n - (N - 1) / 2));
    S_x = SNR * v_a * v_a' + eye(N);
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG_4(k) = abs(w' * v_a)^2 / (w' * w);
    AG_4(k) = 10 * log10(real(AG_4(k)));
end

% Calculate array gain for an SNR = 30 dB

for k = 1 : length(delta)
    SNR = 10^(30 / 10);
    v_a = exp(j * pi * cos(theta) * (1 + delta(k)) * (n - (N - 1) / 2));
    S_x = SNR * v_a * v_a' + eye(N);
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG_5(k) = abs(w' * v_a)^2 / (w' * w);
    AG_5(k) = 10 * log10(real(AG_5(k)));
end

% Plot the results

figure;
plot(delta, AG_1, delta, AG_2, delta, AG_3, delta, AG_4, delta, AG_5);
grid;
xlabel('\Delta = f_a / f_c');
ylabel('Array Gain (dB)');
title('Array Gain: MPDR Beamformer with Frequency Mismatches');
legend('SNR = -10 dB', 'SNR = 0 dB', 'SNR = 10 dB', 'SNR = 20 dB', 'SNR = 30 dB');

