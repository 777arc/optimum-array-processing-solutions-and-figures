% problem6613b.m: Solves Problem 6.6.13(b) in Van Trees, Volume IV.
% 11/17/03, J.A. Tague

% Initialize parameters for a ten-element standard line array

N = 10;
v_s = ones(10, 1);
ind = 0 : N - 1;
n = N - 1 : -1 : 0;
n = n';

B_fI = 0.1;   % The normalized temporal interference bandwidth

% Determine array manifold for the first interference source, assuming that it is narrowband

u_I1 = 0.3;
v1 = exp(j * pi * u_I1 * (n - (N - 1) / 2));

% Generate the actual array covariance matrix due to the first interference source
% at u = 0.3

row1 = exp(j * pi * ind * u_I1) .* sinc(0.5 * ind * B_fI * u_I1);
S_1 = toeplitz(row1);

u_I2 = -1 : 0.005 : 1;

AG1 = zeros(1, length(u_I2));
AG2 = zeros(1, length(u_I2));
AG3 = zeros(1, length(u_I2));
AG4 = zeros(1, length(u_I2));

% Calculate array gain versus u_I2 for INR_2 = -10 dB

INR = -10;

for k = 1 : length(u_I2)
    
    % Calculate the actual array noise covariance matrix using the
    % broadband interference source and uncorrelated noise
    
    row1 = exp(j * pi * ind * u_I2(k)) .* sinc(0.5 * ind * B_fI * u_I2(k));
    S_2 = toeplitz(row1);
    S_ntrue = 100 * S_1 + 10^(INR / 10) * S_2 + eye(10);
    rho_n = N * S_ntrue / trace(S_ntrue);
    
    % Build the optimum MVDR beamformer and calculate array gain, assuming
    % that the noise is narrowband
    
    v2 = exp(j * pi * u_I2(k) * (n - (N - 1) / 2));
    S_ndesign = 100 * v1 * v1' + 10^(INR / 10) * v2 * v2' + eye(10);
    w = inv(S_ndesign) * v_s / (v_s' * inv(S_ndesign) * v_s);
    
    % Calculate the array gain
    
    AG1(k) = abs(w' * v_s)^2 / (w' * rho_n * w);
    
end;

AG1 = 10 * log10(real(AG1));

% Calculate array gain versus u_I2 for INR_2 = 0 dB

INR = 0;

for k = 1 : length(u_I2)
    
    % Calculate the actual array noise covariance matrix using the
    % broadband interference source and uncorrelated noise
    
    row1 = exp(j * pi * ind * u_I2(k)) .* sinc(0.5 * ind * B_fI * u_I2(k));
    S_2 = toeplitz(row1);
    S_ntrue = 100 * S_1 + 10^(INR / 10) * S_2 + eye(10);
    rho_n = N * S_ntrue / trace(S_ntrue);
    
    % Build the optimum MVDR beamformer and calculate array gain, assuming
    % that the noise is narrowband
    
    v2 = exp(j * pi * u_I2(k) * (n - (N - 1) / 2));
    S_ndesign = 100 * v1 * v1' + 10^(INR / 10) * v2 * v2' + eye(10);
    w = inv(S_ndesign) * v_s / (v_s' * inv(S_ndesign) * v_s);
    
    % Calculate the array gain
    
    AG2(k) = abs(w' * v_s)^2 / (w' * rho_n * w);
    
end

AG2 = 10 * log10(real(AG2));

% Calculate array gain versus u_I2 for INR_2 = 10 dB

INR = 10;

for k = 1 : length(u_I2)
    
    % Calculate the actual array noise covariance matrix using the
    % broadband interference source and uncorrelated noise
    
    row1 = exp(j * pi * ind * u_I2(k)) .* sinc(0.5 * ind * B_fI * u_I2(k));
    S_2 = toeplitz(row1);
    S_ntrue = 100 * S_1 + 10^(INR / 10) * S_2 + eye(10);
    rho_n = N * S_ntrue / trace(S_ntrue);
    
    % Build the optimum MVDR beamformer and calculate array gain, assuming
    % that the noise is narrowband
    
    v2 = exp(j * pi * u_I2(k) * (n - (N - 1) / 2));
    S_ndesign = 100 * v1 * v1' + 10^(INR / 10) * v2 * v2' + eye(10);
    w = inv(S_ndesign) * v_s / (v_s' * inv(S_ndesign) * v_s);
    
    % Calculate the array gain
    
    AG3(k) = abs(w' * v_s)^2 / (w' * rho_n * w);
    
end

AG3 = 10 * log10(real(AG3));

% Calculate array gain versus u_I2 for INR_2 = 20 dB

INR = 20;

for k = 1 : length(u_I2)
    
    % Calculate the actual array noise covariance matrix using the
    % broadband interference source and uncorrelated noise
    
    row1 = exp(j * pi * ind * u_I2(k)) .* sinc(0.5 * ind * B_fI * u_I2(k));
    S_2 = toeplitz(row1);
    S_ntrue = 100 * S_1 + 10^(INR / 10) * S_2 + eye(10);
    rho_n = N * S_ntrue / trace(S_ntrue);
    
    % Build the optimum MVDR beamformer and calculate array gain, assuming
    % that the noise is narrowband
    
    v2 = exp(j * pi * u_I2(k) * (n - (N - 1) / 2));
    S_ndesign = 100 * v1 * v1' + 10^(INR / 10) * v2 * v2' + eye(10);
    w = inv(S_ndesign) * v_s / (v_s' * inv(S_ndesign) * v_s);
    
    % Calculate the array gain
    
    AG4(k) = abs(w' * v_s)^2 / (w' * rho_n * w);
    
end

AG4 = 10 * log10(real(AG4));

% Plot the array gain versus u_I2

figure;
plot(u_I2, AG1, u_I2, AG2, u_I2, AG3, u_I2, AG4);
xlabel('Arrival Direction of the Second Interference Source');
ylabel('Array Gain (dB)');
title('Array Gain Versus Arrival Direction of Second Interferer: B_{fI} = 0.1; INR_2 = -10, 0, 10 and 20 dB');
grid;
legend('INR = -10 dB', 'INR = 0 dB', 'INR = 10 dB', 'INR = 20 dB');
    
   

