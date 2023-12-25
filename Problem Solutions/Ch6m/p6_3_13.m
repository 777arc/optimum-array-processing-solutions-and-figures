% Matlab script problem6313.m: Solves Problem 6.3.13 in Van Trees, Volume IV.

% 11/8/03, J.A. Tague

N_x = 5;
U_x = [];
U_y = [];

% Generate the array manifold vector for the signal of interest

vs = ones(19, 1);

% Generate the array manifold vector for the two interference sources

u_R = 0.5;
u_x = u_R * cos(0 * pi / 180);
u_y = u_R * sin(0 * pi / 180);

old_v = [ ];
    
for m = -(N_x - 1) / 2 : (N_x - 1) / 2
    p1 = exp(j * pi * sqrt(3) * m * u_y / 2);
    n = -(N_x - abs(m) - 1)/2 : (N_x - abs(m) - 1)/2;
    n = n';
    p2 = exp(j * pi * n * u_x);
    v = p1 * p2;
    v_man = [v; old_v];
    old_v = v_man;
end

v1 = v_man;

u_R = 0.5;
u_x = u_R * cos(30 * pi / 180);
u_y = u_R * sin(30 * pi / 180);

old_v = [ ];
    
for m = -(N_x - 1) / 2 : (N_x - 1) / 2
    p1 = exp(j * pi * sqrt(3) * m * u_y / 2);
    n = -(N_x - abs(m) - 1)/2 : (N_x - abs(m) - 1)/2;
    n = n';
    p2 = exp(j * pi * n * u_x);
    v = p1 * p2;
    v_man = [v; old_v];
    old_v = v_man;
end

v2 = v_man;

% Calculate S_x assuming the signal power is unity

S_x = vs * vs' + 100 * v1 * v1' + 100 * v2 * v2' + eye(19);

% Calculate the optimum minimum power distortionless response beamformer
% with v_m = vs

w = inv(S_x) * vs / (vs' * inv(S_x) * vs);

% Calculate the array gain

S_n = 100 * v1 * v1' + 100 * v2 * v2' + eye(19);

rho = S_n / (201 + 1);

AG = abs(w' * vs)^2 / (w' * rho * w);
AG = 10 * log10(real(AG))

% Generate and plot the beampattern function

u_x = -1 : 0.01 : 1;
u_y = -1 : 0.01 : 1;

for k1 = 1 : length(u_y)
    for k2 = 1 : length(u_x)
        old_v = [ ];
        for m = -(N_x - 1) / 2 : (N_x - 1) / 2
            p1 = exp(j * pi * sqrt(3) * m * u_y(k1) / 2);
            n = -(N_x - abs(m) - 1)/2 : (N_x - abs(m) - 1)/2;
            n = n';
            p2 = exp(j * pi * n * u_x(k2));
            v = p1 * p2;
            v_man = [v; old_v];
            old_v = v_man;        
        end
        B(k1, k2) = w' * v_man;
    end
end

B = 20 * log10(abs(B));

for k1 = 1 : size(B, 1)
    for k2 = 1 : size(B, 2)
        if (B(k1, k2) < -60) B(k1, k2) = -60;
        end;
    end;
end;

% Plot three dimensional and contour plots of the beam pattern

figure;
[u_x, u_y] = meshgrid(-1 : 0.01 : 1);
mesh(u_x, u_y, B);
grid on;
axis([-1 1 -1 1 -60 0]);
xlabel('u_x', 'Fontsize', 12);
ylabel('u_y', 'Fontsize', 12);
zlabel('Beam Patttern (dB)', 'Fontsize', 12);
title('Beam Pattern of a 19-Element Hexagonal Array with MPDR Beamforming');

figure;
contour(u_x, u_y, B);
xlabel('u_x', 'Fontsize', 12);
ylabel('u_y', 'Fontsize', 12);
title('Beam Pattern Contour Plot: 19-Element Hexagonal Array; MPDR Beamformer');
grid;

hold on;
plot(0.5 * cos(0), 0.5 * sin(0), 'o', 0.5 * cos(30 * pi / 180), 0.5 * sin(30 * pi / 180), 'o');
hold off;

% Calculate and plot the beampattern taken through slices of constant angle

u_R = -1 : 0.005 : 1;

phi = 0 * pi / 180;

for k1 = 1 : length(u_R)
    old_v = [ ];
    for m = -(N_x - 1) / 2 : (N_x - 1) / 2
        p1 = exp(j * pi * sqrt(3) * m * u_R(k1) * sin(phi) / 2);
        n = -(N_x - abs(m) - 1)/2 : (N_x - abs(m) - 1)/2;
        n = n';
        p2 = exp(j * pi * n * u_R(k1) * cos(phi));
        v = p1 * p2;
        v_man = [v; old_v];
        old_v = v_man;        
    end
    B0(k1) = w' * v_man;
end

B0 = 20 * log10(abs(B0));

for k1 = 1 : size(B0, 1)
    if (B0(k1) < -60) B0(k1) = -60;
    end;
end;

figure;
plot(u_R, B0);
axis([-1 1 -60 0]);
xlabel('u_R');
ylabel('Beam Pattern Function (dB');
title('Slice of Beam Pattern Function: \phi = 0 Degrees');
grid;

phi = 10 * pi / 180;

for k1 = 1 : length(u_R)
    old_v = [ ];
    for m = -(N_x - 1) / 2 : (N_x - 1) / 2
        p1 = exp(j * pi * sqrt(3) * m * u_R(k1) * sin(phi) / 2);
        n = -(N_x - abs(m) - 1)/2 : (N_x - abs(m) - 1)/2;
        n = n';
        p2 = exp(j * pi * n * u_R(k1) * cos(phi));
        v = p1 * p2;
        v_man = [v; old_v];
        old_v = v_man;        
    end
    B0(k1) = w' * v_man;
end

B0 = 20 * log10(abs(B0));

for k1 = 1 : size(B0, 1)
    if (B0(k1) < -60) B0(k1) = -60;
    end;
end;

figure;
plot(u_R, B0);
axis([-1 1 -60 0]);
xlabel('u_R');
ylabel('Beam Pattern Function (dB');
title('Slice of Beam Pattern Function: \phi = 10 Degrees');
grid;

phi = 20 * pi / 180;

for k1 = 1 : length(u_R)
    old_v = [ ];
    for m = -(N_x - 1) / 2 : (N_x - 1) / 2
        p1 = exp(j * pi * sqrt(3) * m * u_R(k1) * sin(phi) / 2);
        n = -(N_x - abs(m) - 1)/2 : (N_x - abs(m) - 1)/2;
        n = n';
        p2 = exp(j * pi * n * u_R(k1) * cos(phi));
        v = p1 * p2;
        v_man = [v; old_v];
        old_v = v_man;        
    end
    B0(k1) = w' * v_man;
end

B0 = 20 * log10(abs(B0));

for k1 = 1 : size(B0, 1)
    if (B0(k1) < -60) B0(k1) = -60;
    end;
end;

figure;
plot(u_R, B0);
axis([-1 1 -60 0]);
xlabel('u_R');
ylabel('Beam Pattern Function (dB');
title('Slice of Beam Pattern Function: \phi = 20 Degrees');
grid;

phi = 30 * pi / 180;

for k1 = 1 : length(u_R)
    old_v = [ ];
    for m = -(N_x - 1) / 2 : (N_x - 1) / 2
        p1 = exp(j * pi * sqrt(3) * m * u_R(k1) * sin(phi) / 2);
        n = -(N_x - abs(m) - 1)/2 : (N_x - abs(m) - 1)/2;
        n = n';
        p2 = exp(j * pi * n * u_R(k1) * cos(phi));
        v = p1 * p2;
        v_man = [v; old_v];
        old_v = v_man;        
    end
    B0(k1) = w' * v_man;
end

B0 = 20 * log10(abs(B0));

for k1 = 1 : size(B0, 1)
    if (B0(k1) < -60) B0(k1) = -60;
    end;
end;

figure;
plot(u_R, B0);
axis([-1 1 -60 0]);
xlabel('u_R');
ylabel('Beam Pattern Function (dB');
title('Slice of Beam Pattern Function: \phi = 30 Degrees');
grid;