% problem679b.m: Solves Problem 6.7.9 part (b), Van Trees, Volume IV
% 11/28/03, J.A. Tague

% Initialize parameters for a twenty element standard line array

N = 20;
n = N - 1 : -1 : 0;
n = n';

% Calculate the array manifold vector of the desired signal and interferers

v_m = ones(N, 1);
v_1 = exp(j * pi * 0.3 * (n - (N - 1) / 2));
v_2 = exp(j * pi * 0.5 * (n - (N - 1) / 2));
v_3 = exp(j * pi * 0.7 * (n - (N - 1) / 2));
v_4 = exp(j * pi * -0.5 * (n - (N - 1) / 2));
v_5 = exp(j * pi * -0.3 * (n - (N - 1) / 2));

% Calculate the powers of the first four interferers

sigma_21 = 10^(40 / 10);
sigma_22 = 10^(30 / 10);
sigma_23 = 10^(50 / 10);
sigma_24 = 10^(10 / 10);

% Calculate S_x, the covariance matrix of the four inteferers with known power plus
% white noise

S_n4 = sigma_21 * v_1 * v_1' + sigma_22 * v_2 * v_2' + sigma_23 * v_3 * v_3' + sigma_24 * v_4 * v_4' + eye(N);

% Build the contraint matrix C and gain constraint vector g, including a
% distortionless response constraint at u = 0

u_dir = [0 0.3 0.5 0.7 -0.5];
C = exp(j * pi * (n - (N - 1) / 2) * u_dir);
g = [1 0 0 0 0]';

% Evaluate the performance of the LCMP beamformer implemented as a
% generalized sidelobe canceller

% Calculate the upper branch of the structure and the blocking matrix

w_q = C * inv(C' * C) * g;
[U, Delta] = svd(C);
B = U(:, 6:N);

% Calculate array gain of the MPLC beamformer as a function of the fifth
% interferer's INR.  Assume the DOA's of the four interferers infused into C are 
% known exactly.  This result bounds the processor's performance.

INR_5 = 0 : 5 : 60;
AG_true = zeros(1, length(INR_5));

for k = 1 : length(INR_5)
    sigma_25 = 10^(INR_5(k) / 10);
    S_n = sigma_25 * v_5 * v_5' + S_n4;
    rho_n = (N / trace(S_n)) * S_n;
    S_x = v_m * v_m' + S_n;
    w_a = inv(B' * S_x * B) * B' * S_x' * w_q;
    w = w_q - B * w_a;
    AG = abs(w' * v_m)^2 / (w' * rho_n * w);
    AG_true(k) = 10 * log10(real(AG));
end

% Estimate the expected value of array gain assuming random errors in the
% interferer's estimated arrival directions

N_iter = 1000;

sigma_e = 0.1;
AG_est1 = zeros(1, length(INR_5));

for k = 1 : length(INR_5)
    
    AG_hat = 0;
    sigma_25 = 10^(INR_5(k) / 10);
    
    for m = 1 : N_iter
        u_dir = [0 0.3 0.5 0.7 -0.5] + [0 sigma_e * randn(1, 4)];
        C = exp(j * pi * (n - (N - 1) / 2) * u_dir);
        w_q = C * inv(C' * C) * g;
        [U, Delta] = svd(C);
        B = U(:, 6:N);
        S_n = sigma_25 * v_5 * v_5' + S_n4;
        rho_n = (N / trace(S_n)) * S_n;
        S_x = v_m * v_m' + S_n;
        w_a = inv(B' * S_x * B) * B' * S_x' * w_q;
        w = w_q - B * w_a;
        AG_hat = abs(w' * v_m)^2 / (w' * rho_n * w) + AG_hat;
    end;
    
    AG_hat = AG_hat / N_iter;
    AG_est1(k) = 10 * log10(real(AG_hat));
    
end;

sigma_e = 0.2;
AG_est2 = zeros(1, length(INR_5));

for k = 1 : length(INR_5)
    
    AG_hat = 0;
    sigma_25 = 10^(INR_5(k) / 10);
    
    for m = 1 : N_iter
        u_dir = [0 0.3 0.5 0.7 -0.5] + [0 sigma_e * randn(1, 4)];
        C = exp(j * pi * (n - (N - 1) / 2) * u_dir);
        w_q = C * inv(C' * C) * g;
        [U, Delta] = svd(C);
        B = U(:, 6:N);
        S_n = sigma_25 * v_5 * v_5' + S_n4;
        rho_n = (N / trace(S_n)) * S_n;
        S_x = v_m * v_m' + S_n;
        w_a = inv(B' * S_x * B) * B' * S_x' * w_q;
        w = w_q - B * w_a;
        AG_hat = abs(w' * v_m)^2 / (w' * rho_n * w) + AG_hat;
    end;
    
    AG_hat = AG_hat / N_iter;
    AG_est2(k) = 10 * log10(real(AG_hat));
    
end;

sigma_e = 0.3;
AG_est3 = zeros(1, length(INR_5));

for k = 1 : length(INR_5)
    
    AG_hat = 0;
    sigma_25 = 10^(INR_5(k) / 10);
    
    for m = 1 : N_iter
        u_dir = [0 0.3 0.5 0.7 -0.5] + [0 sigma_e * randn(1, 4)];
        C = exp(j * pi * (n - (N - 1) / 2) * u_dir);
        w_q = C * inv(C' * C) * g;
        [U, Delta] = svd(C);
        B = U(:, 6:N);
        S_n = sigma_25 * v_5 * v_5' + S_n4;
        rho_n = (N / trace(S_n)) * S_n;
        S_x = v_m * v_m' + S_n;
        w_a = inv(B' * S_x * B) * B' * S_x' * w_q;
        w = w_q - B * w_a;
        AG_hat = abs(w' * v_m)^2 / (w' * rho_n * w) + AG_hat;
    end;
    
    AG_hat = AG_hat / N_iter;
    AG_est3(k) = 10 * log10(real(AG_hat));
    
end;

% Plot the results

floor = 0;

for k = 1 : length(INR_5)
    if (AG_true(k) < floor) AG_true(k) = floor;
    end;
    if (AG_est1(k) < floor) AG_est1(k) = floor;
    end;
    if (AG_est2(k) < floor) AG_est2(k) = floor;
    end;
    if (AG_est3(k) < floor) AG_est3(k) = floor;
    end;
end;

figure;
plot(INR_5, AG_true, INR_5, AG_est1, INR_5, AG_est2, INR_5, AG_est3);
grid;
xlabel('INR of the Fifth Interferer at u_I = -0.3 (dB)');
ylabel('Array Gain (dB)');
title('Array Gain of the LCMP-GSC Beamformer: Known and Estimated Interferer Arrival Directions');
legend('Known Interferer DOA', 'Estimated Interferer DOA, \sigma_e = 0.1', 'Estimated Interferer DOA, \sigma_e = 0.2', 'Estimated Interferer DOA: \sigma_e = 0.3');