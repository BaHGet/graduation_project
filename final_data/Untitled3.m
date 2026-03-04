%% Hexagonal Lattice Band Structure Calculator (EPWE Method)
% Based on the Empty Plane Wave Expansion (EPWE) approach
% Calculates and visualizes band structure for a 2D hexagonal metallic lattice
%
% Steps:
%   1. Construct hexagonal geometry with potential barriers
%   2. Define reciprocal lattice vectors
%   3. Build and solve the Hamiltonian (EPWE)
%   4. Plot the lattice and band structure along high-symmetry path

clear; clc; close all;

%% =========================================================
%  SECTION 1 — PROBLEM INPUTS
%  Material & lattice parameters
%% =========================================================

% Lattice constant (nm)
a = 1.0;

% Effective mass ratio (m*/m_e) — metallic-like material
m_eff_ratio = 0.4;

% Physical constants
hbar   = 1.0546e-34;   % J·s
m_e    = 9.109e-31;    % kg  (free electron mass)
eV2J   = 1.602e-19;    % 1 eV in Joules
m_eff  = m_eff_ratio * m_e;

% Potential parameters
V0     = 0.5;          % Barrier height (eV)  — potential inside the discs
r_disc = 0.30 * a;     % Disc radius (fraction of lattice constant)

% Number of plane waves (reciprocal lattice vectors) to include
% Nmax sets the shell index; total PWs ? 3*Nmax^2 + 3*Nmax + 1
Nmax   = 3;

% Number of bands to display
n_bands = 6;

%% =========================================================
%  SECTION 2 — CONSTRUCT HEXAGONAL GEOMETRY
%% =========================================================

% --- Primitive lattice vectors (2D hexagonal) ---
a1 = a * [1,        0      ];
a2 = a * [1/2, sqrt(3)/2   ];

% --- Reciprocal lattice vectors ---
% For 2D: b1, b2 satisfy ai · bj = 2? ?_ij
Area = abs(a1(1)*a2(2) - a1(2)*a2(1));   % unit cell area
b1   = (2*pi/Area) * [ a2(2), -a2(1)];
b2   = (2*pi/Area) * [-a1(2),  a1(1)];

fprintf('Lattice constant a = %.2f nm\n', a);
fprintf('Reciprocal vectors: |b1| = %.4f nm^-1, |b2| = %.4f nm^-1\n', ...
        norm(b1), norm(b2));

% --- Generate set of reciprocal lattice vectors G ---
G_list = [];
for n1 = -Nmax:Nmax
    for n2 = -Nmax:Nmax
        G = n1*b1 + n2*b2;
        G_list = [G_list; G, n1, n2]; %#ok<AGROW>
    end
end
% Keep only vectors within a reasonable cutoff
G_mag   = sqrt(G_list(:,1).^2 + G_list(:,2).^2);
G_cut   = (Nmax + 0.5) * norm(b1);
keep    = G_mag <= G_cut;
G_list  = G_list(keep, :);
N_G     = size(G_list, 1);
fprintf('Number of plane waves used: %d\n', N_G);

%% =========================================================
%  SECTION 3 — POTENTIAL FORM FACTORS  V(G-G')
%  Circular disc potential: V(r) = V0 inside radius r_disc
%  Fourier component: V(G) = V0 * f * 2*J1(|G|*r_disc)/(|G|*r_disc)
%  where f = filling fraction = ? r_disc^2 / Area
%% =========================================================

f_fill = pi * r_disc^2 / Area;   % filling fraction

% Pre-compute all needed V(G-G') — store in matrix
V_mat = zeros(N_G, N_G);
for i = 1:N_G
    for j = 1:N_G
        dG = G_list(i,1:2) - G_list(j,1:2);
        q  = norm(dG);
        if q < 1e-10
            % G=G' term
            V_mat(i,j) = V0 * f_fill;
        else
            % Cylindrical form factor
            V_mat(i,j) = V0 * f_fill * 2 * besselj(1, q*r_disc) / (q*r_disc);
        end
    end
end
% Convert to Joules for Hamiltonian
V_mat_J = V_mat * eV2J;

%% =========================================================
%  SECTION 4 — HIGH-SYMMETRY k-PATH  (? ? M ? K ? ?)
%  Standard hexagonal BZ path
%% =========================================================

Gamma = [0, 0];
M     = 0.5 * b1;
K     = (b1 + b2) / 3 * 2;   % K point at 2/3 (b1+b2)/2... recalculate:
K     = (1/3)*b1 + (1/3)*b2 + (1/3)*(b1-b2);
% Correct high-symmetry points for hexagonal lattice
% ?=(0,0), M=b1/2, K=(b1+b2)/3 * ... use standard definition
M_pt  = (b1 + b2) / 2;       % M point (midpoint of BZ edge)
K_pt  = (2*b1 + b2) / 3;     % K point (BZ corner)

n_seg = 50;  % k-points per segment

k_GM  = kpath_segment(Gamma, M_pt, n_seg);
k_MK  = kpath_segment(M_pt,  K_pt,  n_seg);
k_KG  = kpath_segment(K_pt,  Gamma, n_seg);

k_path    = [k_GM; k_MK(2:end,:); k_KG(2:end,:)];
N_k       = size(k_path, 1);

% Build x-axis (cumulative distance)
k_dist = zeros(N_k, 1);
for ik = 2:N_k
    k_dist(ik) = k_dist(ik-1) + norm(k_path(ik,:) - k_path(ik-1,:));
end

% Tick positions for high-symmetry labels
tick_pos  = [0, k_dist(n_seg), k_dist(2*n_seg-1), k_dist(end)];

%% =========================================================
%  SECTION 5 — EPWE: BUILD & DIAGONALISE HAMILTONIAN
%% =========================================================

% Pre-factor for kinetic energy: ?²/(2 m_eff) in eV·nm²
KE_pref = (hbar^2 / (2 * m_eff)) / eV2J * (1e9)^2;  % eV·nm²
% (G vectors are in nm^-1 if a in nm)

E_bands = zeros(N_k, N_G);

fprintf('Solving eigenvalue problem for %d k-points...\n', N_k);
for ik = 1:N_k
    kv = k_path(ik, :);   % current k-vector (nm^-1)

    % Kinetic energy diagonal: ?²|k+G|²/(2m*)
    T_diag = zeros(N_G, 1);
    for ig = 1:N_G
        kG     = kv + G_list(ig, 1:2);
        T_diag(ig) = KE_pref * (kG(1)^2 + kG(2)^2);
    end

    % Full Hamiltonian H = T + V
    H = V_mat + diag(T_diag);

    % Solve (H is Hermitian)
    evals = sort(real(eig((H + H') / 2)));
    E_bands(ik, :) = evals;
end
fprintf('Done.\n');

%% =========================================================
%  SECTION 6 — FIGURE 1: HEXAGONAL LATTICE WITH POTENTIAL
%% =========================================================

figure('Name','Hexagonal Lattice — Potential Landscape', ...
       'Color','w','Position',[50 100 560 500]);

% Draw superlattice: 5×5 unit cells
n_cells = 3;
hold on; axis equal; box on;

% Background colour
ax = gca;
ax.Color = [0.95 0.97 1.0];
set(ax,'XTick',[],'YTick',[]);

% Draw potential discs (barrier regions)
theta_circ = linspace(0, 2*pi, 120);
xc = r_disc * cos(theta_circ);
yc = r_disc * sin(theta_circ);

disc_color  = [0.15 0.15 0.55];   % dark blue for barrier
metal_color = [0.85 0.90 0.98];   % light blue for free region

% Fill background (free-electron region)
xlim_val = n_cells * a * 1.8;
fill([-xlim_val xlim_val xlim_val -xlim_val], ...
     [-xlim_val -xlim_val xlim_val xlim_val], metal_color, ...
     'EdgeColor','none');

% Place discs on every hexagonal lattice site
for n1 = -n_cells:n_cells
    for n2 = -n_cells:n_cells
        site = n1*a1 + n2*a2;
        % Rough visibility filter
        if abs(site(1)) > xlim_val*1.1 || abs(site(2)) > xlim_val*1.1
            continue;
        end
        fill(site(1)+xc, site(2)+yc, disc_color, ...
             'EdgeColor',[0 0 0.3],'LineWidth',0.5,'FaceAlpha',0.85);
    end
end

% Draw unit cell outline
cell_verts = [0 0; a1; a1+a2; a2];
patch('XData',cell_verts(:,1),'YData',cell_verts(:,2), ...
      'FaceColor','none','EdgeColor',[0.8 0.2 0.2],'LineWidth',2);

% Draw lattice vectors
quiver(0,0, a1(1),a1(2), 0,'r','LineWidth',2,'MaxHeadSize',0.5);
quiver(0,0, a2(1),a2(2), 0,'g','LineWidth',2,'MaxHeadSize',0.5);
text(a1(1)+0.05, a1(2),       '\bfa_1','Color','r','FontSize',13);
text(a2(1)-0.18, a2(2)+0.05,  '\bfa_2','Color',[0 0.55 0],'FontSize',13);

xlabel('x (nm)','FontSize',12);
ylabel('y (nm)','FontSize',12);
title({'2D Hexagonal Lattice — Potential Landscape'; ...
       sprintf('V_0 = %.2f eV,  r_{disc}/a = %.2f,  a = %.1f nm', ...
               V0, r_disc/a, a)}, 'FontSize',13);

% Colorbar annotation
ax2 = axes('Position',[0.82 0.30 0.03 0.40],'Visible','off');
colormap(ax2, [metal_color; disc_color]);
cb = colorbar(ax2,'Location','east');
cb.Ticks      = [0 1];
cb.TickLabels = {'0', sprintf('%.1f eV',V0)};
cb.Label.String = 'V_s(r)';
cb.Label.FontSize = 11;

hold off;

%% =========================================================
%  SECTION 7 — FIGURE 2: BAND STRUCTURE
%% =========================================================

figure('Name','Band Structure — Hexagonal Lattice (EPWE)', ...
       'Color','w','Position',[650 100 620 520]);

hold on;
colors = lines(n_bands);

for ib = 1:min(n_bands, N_G)
    plot(k_dist, E_bands(:, ib), '-', 'Color', colors(ib,:), 'LineWidth', 2);
end

% High-symmetry vertical lines
for tp = tick_pos
    xline(tp, '--k', 'LineWidth', 0.8, 'Alpha', 0.4);
end

ax = gca;
ax.XTick      = tick_pos;
ax.XTickLabel = {'\Gamma', 'M', 'K', '\Gamma'};
ax.FontSize   = 13;

xlim([0, k_dist(end)]);
ylim_max = max(E_bands(:, min(n_bands, N_G))) * 1.08;
ylim_min = max(0, min(E_bands(:,1)) - 0.05);
ylim([ylim_min, ylim_max]);

xlabel('k-path', 'FontSize',14);
ylabel('Electron Energy (eV)', 'FontSize',14);
title({'Band Structure — 2D Hexagonal Metallic Lattice (EPWE)'; ...
       sprintf('m^* = %.1f m_e,  V_0 = %.2f eV,  a = %.1f nm', ...
               m_eff_ratio, V0, a)}, 'FontSize',13);

legend(arrayfun(@(b) sprintf('Band %d',b), 1:min(n_bands,N_G), ...
       'UniformOutput',false), 'Location','northeast','FontSize',10);

box on; grid on; ax.GridAlpha = 0.2;
hold off;