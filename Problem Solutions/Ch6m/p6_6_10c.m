% problem6610c.m: Solves Problem 6.6.10 part (c) in Van Trees, Volume IV
% 11/16/03, J.A. Tague

n_iter = 20000; % Number of iterations in Monte Carlo trials to estimate E{array gain}

% Initialize parameters for a ten element standard line array

N = 10;
n = N - 1 : -1 : 0;
n = n';

% Compute the array manifold vector for the optimum beamformer

theta = 30 * pi / 180;
v_m = exp(j * pi * cos(theta) * (n - (N - 1) / 2));

SNR_db = -10 : 5 : 30;

AG_1 = zeros(1, length(SNR_db));
AG_2 = zeros(1, length(SNR_db));
AG_3 = zeros(1, length(SNR_db));

% Estimate E{Array Gain} for f_1 = 0.04 * f_c

f_1 = 0.04;

for k = 1 : length(SNR_db)
    
    SNR = 10^(SNR_db(k) / 10);
    AG_hat = 0;
    M = N * SNR;
    
    for m = 1 : n_iter
        delta = 2 * f_1 * (rand(1) - 0.5);
        v_a = exp(j * pi * (1 + delta) * cos(theta) * (n - (N - 1)/2));
        rho = v_m' * v_a / N;
        AG_hat = AG_hat + (N * abs(rho)^2) / (1 + (2 * M + M^2) * (1 - abs(rho)^2));
    end;
    
    AG_hat = AG_hat / n_iter;
    AG_1(k) = 10 * log10(AG_hat);
    
end;

% Estimate E{Array Gain} for f_1 = 0.1 * f_c

f_1 = 0.1;

for k = 1 : length(SNR_db)
    
    SNR = 10^(SNR_db(k) / 10);
    AG_hat = 0;
    M = N * SNR;
    
    for m = 1 : n_iter
        delta = 2 * f_1 * (rand(1) - 0.5);
        v_a = exp(j * pi * (1 + delta) * cos(theta) * (n - (N - 1)/2));
        rho = v_m' * v_a / N;
        AG_hat = AG_hat + (N * abs(rho)^2) / (1 + (2 * M + M^2) * (1 - abs(rho)^2));
    end;
    
    AG_hat = AG_hat / n_iter;
    AG_2(k) = 10 * log10(AG_hat);
    
end;

% Estimate E{Array Gain} for f_1 = 0.2 * f_c

f_1 = 0.2;

for k = 1 : length(SNR_db)
    
    SNR = 10^(SNR_db(k) / 10);
    AG_hat = 0;
    M = N * SNR;
    
    for m = 1 : n_iter
        delta = 2 * f_1 * (rand(1) - 0.5);
        v_a = exp(j * pi * (1 + delta) * cos(theta) * (n - (N - 1)/2));
        rho = v_m' * v_a / N;
        AG_hat = AG_hat + (N * abs(rho)^2) / (1 + (2 * M + M^2) * (1 - abs(rho)^2));
    end;
    
    AG_hat = AG_hat / n_iter;
    AG_3(k) = 10 * log10(AG_hat);
    
end;

% Plot the results

figure;
plot(SNR_db, AG_1, SNR_db, AG_2, SNR_db, AG_3);
grid;
xlabel('Signal-to-Noise Ratio (dB)');
ylabel('Estimated E\{Array Gain\}');
title('Estimated Expected Value of Array Gain: MPDR Beamformer');
legend('f_1 = 0.04 f_c', 'f_1 = 0.1 f_c', 'f_1 = 0.2 f_c');

