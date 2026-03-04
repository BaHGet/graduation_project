clc;
clear;
close all;

%% =========================
% PARAMETERS
%% =========================

a = 1;              % lattice constant (normalized)
R = 0.25;           % radius of potential barrier
Nx = 4;             % number of cells in x
Ny = 4;             % number of cells in y

hbar = 1;
m_eff = 1;
V0 = 5;             % barrier height

%% =========================
% HEXAGONAL LATTICE VECTORS
%% =========================

a1 = a * [1 0];
a2 = a * [1/2 sqrt(3)/2];

%% =========================
% DRAW LATTICE
%% =========================

figure;
hold on;
axis equal;
title('Hexagonal Lattice with Circular Potential Barriers');

theta = linspace(0,2*pi,100);

for n1 = -Nx:Nx
    for n2 = -Ny:Ny
        
        center = n1*a1 + n2*a2;
        
        x = center(1) + R*cos(theta);
        y = center(2) + R*sin(theta);
        
        fill(x,y,[0.8 0.8 0.8],'EdgeColor','k');
        
    end
end

xlabel('x');
ylabel('y');
grid on;

%% =========================
% RECIPROCAL LATTICE
%% =========================

area = a1(1)*a2(2) - a1(2)*a2(1);
b1 = (2*pi/area) * [ a2(2) -a2(1)];
b2 = (2*pi/area) * [-a1(2)  a1(1)];

%% =========================
% PLANE WAVE BASIS
%% =========================

Nmax = 2;
G = [];

for n1 = -Nmax:Nmax
    for n2 = -Nmax:Nmax
        G = [G; n1*b1 + n2*b2];
    end
end

NG = size(G,1);

%% =========================
% HIGH SYMMETRY POINTS
%% =========================

Gamma = [0 0];
K = (b1 + 2*b2)/3;
M = (b1 + b2)/2;

Nk_segment = 40;

k_path = [];

% Gamma -> K
for i=0:Nk_segment-1
    t=i/(Nk_segment-1);
    k_path=[k_path; (1-t)*Gamma + t*K];
end

% K -> M
for i=0:Nk_segment-1
    t=i/(Nk_segment-1);
    k_path=[k_path; (1-t)*K + t*M];
end

% M -> Gamma
for i=0:Nk_segment-1
    t=i/(Nk_segment-1);
    k_path=[k_path; (1-t)*M + t*Gamma];
end

Nk = size(k_path,1);

%% =========================
% HAMILTONIAN + BAND STRUCTURE
%% =========================

Ebands = zeros(NG,Nk);

for ik = 1:Nk
    
    k = k_path(ik,:);
    H = zeros(NG);
    
    for i = 1:NG
        Gi = G(i,:);
        H(i,i) = (hbar^2/(2*m_eff))*norm(k+Gi)^2;
    end
    
    % simple diagonal potential
    H = H + V0*eye(NG);
    
    E = eig(H);
    Ebands(:,ik) = sort(real(E));
end

%% =========================
% PLOT BAND STRUCTURE
%% =========================

figure;
hold on;

for n = 1:6
    plot(Ebands(n,:),'LineWidth',1.5);
end

set(gca,'XTick',[1 Nk_segment 2*Nk_segment 3*Nk_segment]);
set(gca,'XTickLabel',{'?','K','M','?'});

ylabel('Energy');
title('Band Structure - Hexagonal Lattice');
grid on;