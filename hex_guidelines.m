clear; clf;

% === Parameters ===
a = 1;            % length of vector a
b = 1;            % length of vector b
phi = 90;        % angle between a and b (degrees)
phi_rad = phi * (pi / 180);   % convert degrees ? radians

% === Create hexagon coordinates ===
theta = (0:60:360) * pi/180;
R = 1.0;  % hexagon radius
hex_x = R * cos(theta);
hex_y = R * sin(theta);

% === Plot hexagon ===
plot(hex_x, hex_y, 'b--', 'LineWidth', 1.5); hold on;
axis equal; axis off;

% === Find bottom-left corner position ===
corner_x = R * cos(4*pi/3);
corner_y = R * sin(4*pi/3);

% === Define lattice vectors starting from bottom-left corner ===
origin = [corner_x, corner_y];
va = [a, 0];  % along x
vb = [b*cos(phi_rad), b*sin(phi_rad)];

% === Plot the two vectors ===
quiver(origin(1), origin(2), va(1), va(2), 0, ...
    'Color',[0.7 0.1 0.4],'LineWidth',3,'MaxHeadSize',0.4);
quiver(origin(1), origin(2), vb(1), vb(2), 0, ...
    'Color',[0.7 0.1 0.4],'LineWidth',3,'MaxHeadSize',0.4);

% === Add vector labels ===
text(origin(1)+va(1)*1.05, origin(2)+va(2)*1.05, 'a', ...
     'Color',[0.7 0.1 0.9], 'FontSize',14,'FontWeight','bold');
text(origin(1)+vb(1)*0.89, origin(2)+vb(2)*1.05, 'b', ...
     'Color',[0.7 0.1 0.9], 'FontSize',14,'FontWeight','bold');

% === Draw the ? angle arc ===
r = 0.25 * a; % radius for the angle arc
ang = linspace(0, phi_rad, 100);
arc_x = origin(1) + r*cos(ang);
arc_y = origin(2) + r*sin(ang);
plot(arc_x, arc_y, 'r', 'LineWidth', 1.2);

% === Add ? label ===
text(origin(1)+0.35*a*cos(phi_rad/2), ...
     origin(2)+0.35*a*sin(phi_rad/2)+0.05, ...
     '\phi', 'Color',[0.7 0.1 0.9], 'FontSize',14,'FontWeight','bold');

title('Lattice vectors at the bottom-left corner of hexagon','FontSize',12);
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
