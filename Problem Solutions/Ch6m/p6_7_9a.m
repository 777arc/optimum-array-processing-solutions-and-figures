% problem679a.m: Solves Problem 6.7.9 part (a), Van Trees, Volume IV
% 11/27/03, J.A. Tague

% Initialize parameters for a twenty element standard line array

N = 20;
n = N - 1 : -1 : 0;
n = n';

% Calculate the array manifold vector of the desired signal and the five
% interferers

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

% Calculate S_x, the covariance matrix of the four known inteferers plus
% white noise

S_n4 = sigma_21 * v_1 * v_1' + sigma_22 * v_2 * v_2' + sigma_23 * v_3 * v_3' + sigma_24 * v_4 * v_4' + eye(N);

% Build the contraint matrix C and gain constraint vector g, including a
% distortionless response constraint at u = 0

u_dir = [0 0.3 0.5 0.7 -0.5];
C = exp(j * pi * (n - (N - 1) / 2) * u_dir);
g = [1 0 0 0 0]';

% Calculate array gain of the MPLC beamformer as a function of the fifth
% interferer's INR

INR_5 = 0 : 60;

AGLCMP = zeros(1, length(INR_5));
AGMPDR = zeros(1, length(INR_5));
AGGSC  = zeros(1, length(INR_5));

% Evaluate the performance of the MPDR beamformer

for k = 1 : length(INR_5)
    sigma_25 = 10^(INR_5(k) / 10);
    S_n = sigma_25 * v_5 * v_5' + S_n4;
    rho_n = (N / trace(S_n)) * S_n;
    S_x = v_m * v_m' + S_n;
    w = inv(S_x) * v_m / (v_m' * inv(S_x) * v_m);
    AG = abs(w' * v_m)^2 / (w' * rho_n * w);
    AGMPDR(k) = 10 * log10(real(AG));
end

% Evaluate the performance of the LCMP beamformer implmented in its direct
% form

for k = 1 : length(INR_5)
    sigma_25 = 10^(INR_5(k) / 10);
    S_n = sigma_25 * v_5 * v_5' + S_n4;
    rho_n = (N / trace(S_n)) * S_n;
    S_x = v_m * v_m' + S_n;
    w = inv(S_x) * C * inv(C' * inv(S_x) * C) * g;
    AG = abs(w' * v_m)^2 / (w' * rho_n * w);
    AGLCMP(k) = 10 * log10(real(AG));
end

% Evaluate the performance of the LCMP beamformer implemented as a
% generalized sidelobe canceller

% Calculate the upper branch of the structure

w_q = C * inv(C' * C) * g;

% Calculate the blocking matrix

[U, Delta] = svd(C);
B = U(:, 6:N);

for k = 1 : length(INR_5)
    sigma_25 = 10^(INR_5(k) / 10);
    S_n = sigma_25 * v_5 * v_5' + S_n4;
    rho_n = (N / trace(S_n)) * S_n;
    S_x = v_m * v_m' + S_n;
    w_a = inv(B' * S_x * B) * B' * S_x' * w_q;
    w = w_q - B * w_a;
    AG = abs(w' * v_m)^2 / (w' * rho_n * w);
    AGGSC(k) = 10 * log10(real(AG));
end

% Plot the results

floor = 0;

for k = 1 : length(INR_5)
    if (AGLCMP(k) < floor) AGLCMP(k) = floor;
    end;
    if (AGMPDR(k) < floor) AGMPDR(k) = floor;
    end;
    if (AGGSC(k) < floor)  AGGSC(k) = floor;
    end;
end;

figure;
plot(INR_5, AGLCMP, INR_5, AGMPDR, INR_5, AGGSC);
grid;
xlabel('INR of the Fifth Interferer at u_I = -0.3 (dB)');
ylabel('Array Gain (dB)');
title('Array Gain of the MPLC and MPDR Beamformers');
legend('LCMP Beamformer: Direct Form', 'MPDR Beamformer', 'LCMP Beamformer: GSC');