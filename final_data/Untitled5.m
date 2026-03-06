% =========================================================================
% EPWE Calculation for Hexagonal Lattice Band Structure
% =========================================================================
clear; clc; close all;

%% 1. Problem Inputs & Itemize Region Properties
% These are representative parameters. You can adjust them to match your 
% specific metallic material's effective mass and potential landscape.
a = 1.0;            % Lattice constant (e.g., in nm)
R = 0.3 * a;        % Radius of the circular potential regions
V0 = -2.5;          % Potential depth/barrier (eV)
m_eff = 1.0;        % Effective mass (relative to free electron mass m0)
C_kin = 3.80998 / m_eff; % hbar^2 / 2m0 scaling factor in eV*nm^2

%% 2. Construct Geometry (Real Space - Hexagonal Shape)
% Hexagonal real space translation vectors
a1 = a * [1, 0];
a2 = a * [0.5, sqrt(3)/2];

figure('Name', 'Hexagonal Cluster Geometry'); hold on;
N_rings = 3; % Number of rings to form the large hexagon

for n1 = -N_rings:N_rings
    % Calculate the range for n2 to maintain a hexagonal boundary
    % Condition: |n1| <= N, |n2| <= N, and |n1 + n2| <= N
    n2_min = max(-N_rings, -n1 - N_rings);
    n2_max = min(N_rings, -n1 + N_rings);
    
    for n2 = n2_min:n2_max
        % Calculate real space lattice point center
        R_vec = n1 * a1 + n2 * a2;
        
        % Draw circular potential barrier
        theta = linspace(0, 2*pi, 50);
        x_circle = R_vec(1) + R * cos(theta);
        y_circle = R_vec(2) + R * sin(theta);
        
        % Plotting the circle outline
        plot(x_circle, y_circle, 'k-', 'LineWidth', 1.2);
        
        % Fill the circle to represent the potential Vs(r) regions
        fill(x_circle, y_circle, [0.8 0.8 0.8], 'EdgeColor', 'none'); 
    end
end

axis equal; 
grid off; 
title(['Constructed Hexagonal Cluster (N = ', num2str(N_rings), ')']);
xlabel('x (nm)'); ylabel('y (nm)');
set(gca, 'Color', 'w'); % Set background to white

%% 3. Numerical Considerations (Reciprocal Space & Discretization)
% Reciprocal lattice vectors for hexagonal lattice
b1 = (2*pi/a) * [1, -1/sqrt(3)];
b2 = (2*pi/a) * [0, 2/sqrt(3)];

% Optimize discretization: Define finite reciprocal lattice vectors (G)
N_G_max = 3; % Cutoff for G-vectors (increase for higher accuracy, 3 yields ~49 vectors)
G_vectors = [];
for n1 = -N_G_max:N_G_max
    for n2 = -N_G_max:N_G_max
        G_vectors = [G_vectors; n1*b1 + n2*b2];
    end
end
num_G = size(G_vectors, 1);

%% 4. Define Potential V(G) Matrix Elements
% Calculate filling fraction of the circular region
A_cell = a^2 * sqrt(3) / 2; % Area of hexagonal unit cell
f = pi * R^2 / A_cell;

% Create Fourier coefficient matrix for the potential V(G - G')
V_mat = zeros(num_G, num_G);
for i = 1:num_G
    for j = 1:num_G
        G_diff = G_vectors(i, :) - G_vectors(j, :);
        G_mag = norm(G_diff);
        
        % Analytical Fourier transform of a 2D circular step potential
        if G_mag < 1e-6 
            V_mat(i, j) = V0 * f; % G = 0 case
        else
            V_mat(i, j) = V0 * f * 2 * besselj(1, G_mag*R) / (G_mag*R);
        end
    end
end

%% 5. Define k-point Path (Brillouin Zone High Symmetry Points)
% Standard high-symmetry points for a hexagonal lattice
Gamma = [0, 0];
M     = (pi/a) * [1, 1/sqrt(3)];
K     = (4*pi/(3*a)) * [1, 0];

num_k_pts = 50; % Resolution of the band structure plot
k1 = [linspace(Gamma(1), M(1), num_k_pts)', linspace(Gamma(2), M(2), num_k_pts)'];
k2 = [linspace(M(1), K(1), num_k_pts)',     linspace(M(2), K(2), num_k_pts)'];
k3 = [linspace(K(1), Gamma(1), num_k_pts)', linspace(K(2), Gamma(2), num_k_pts)'];

% Combine path (removing duplicate joining points)
k_path = [k1(1:end-1, :); k2(1:end-1, :); k3];
total_k = size(k_path, 1);

%% 6. EPWE Calculations: Solve the Eigenvalue Problem
num_bands = 4; % Number of lowest energy bands to extract
E_bands = zeros(total_k, num_bands);

% Loop over all k-points in the path
for i = 1:total_k
    k = k_path(i, :);
    H = zeros(num_G, num_G);
    
    % Construct the Hamiltonian matrix
    for row = 1:num_G
        for col = 1:num_G
            if row == col
                % Diagonal elements: Kinetic Energy + V(0)
                k_plus_G = k + G_vectors(row, :);
                H(row, col) = C_kin * norm(k_plus_G)^2 + V_mat(row, col);
            else
                % Off-diagonal elements: V(G - G')
                H(row, col) = V_mat(row, col);
            end
        end
    end
    
    % Solve the eigenvalue problem H * c = E * c
    evals = eig(H);
    evals = sort(real(evals)); % Sort energies from lowest to highest
    E_bands(i, :) = evals(1:num_bands);
end

%% 7. Problem Outputs: Plot the Band Structure
figure('Name', 'Band Structure E_nk'); hold on;
x_axis = 1:total_k;

% Plot using black dots to match the style of your third image
for b = 1:num_bands
    plot(x_axis, E_bands(:, b), 'k.', 'MarkerSize', 6);
end

% Annotate high symmetry points on the X-axis
x_M = num_k_pts;
x_K = 2*num_k_pts - 1;
x_Gamma_end = total_k;

% Get current Y-axis limits to draw vertical lines spanning the whole graph
y_limits = ylim; 

% Draw vertical lines for high-symmetry points (Compatible with older MATLAB versions)
plot([1, 1], y_limits, 'k-', 'HandleVisibility', 'off');
plot([x_M, x_M], y_limits, 'k-', 'HandleVisibility', 'off');
plot([x_K, x_K], y_limits, 'k-', 'HandleVisibility', 'off');
plot([x_Gamma_end, x_Gamma_end], y_limits, 'k-', 'HandleVisibility', 'off');

% Set custom X-axis ticks and labels
set(gca, 'XTick', [1, x_M, x_K, x_Gamma_end]);
set(gca, 'XTickLabel', {'\Gamma', 'M', 'K', '\Gamma'});

ylabel('Electron Energy (eV)');
title('Calculated Band Structure');
box on;