function draw_lattice_cell(a, b, phi)
    % draw_lattice_cell(a, b, phi)
    % Draws a parallelogram (unit cell) defined by lattice vectors a and b
    % with angle phi (in degrees).

    clf; hold on; axis equal; axis off;
phi
    % Convert phi to radians
    phi_rad = phi;% * (pi / 180);

    % Define lattice vectors
    va = [a, 0];
    vb = [b*cos(phi_rad), b*sin(phi_rad)];

    % Define the corners of the cell (parallelogram)
    cell_x = [0, va(1), va(1)+vb(1), vb(1), 0];
    cell_y = [0, va(2), va(2)+vb(2), vb(2), 0];

    % Draw the cell
    plot(cell_x, cell_y, 'w--', 'LineWidth', 1.5);

    % Draw the lattice vectors
    quiver(0, 0, va(1), va(2), 0, 'Color',[0.7 0.1 0.4], 'LineWidth', 3, 'MaxHeadSize', 0.4);
    quiver(0, 0, vb(1), vb(2), 0, 'Color',[0.7 0.1 0.4], 'LineWidth', 3, 'MaxHeadSize', 0.4);

    % Add labels
    text(va(1)*1.05, va(2)*1.05, 'a', 'Color',[0.7 0.1 0.9], 'FontSize', 14, 'FontWeight', 'bold');
    text(vb(1)*1.05, vb(2)*1.05, 'b', 'Color',[0.7 0.1 0.9], 'FontSize', 14, 'FontWeight', 'bold');

    % Draw the phi angle arc
    r = 0.3 * a; % radius for the angle arc
    ang = linspace(0, phi_rad, 100);
    arc_x = r * cos(ang);
    arc_y = r * sin(ang);
    plot(arc_x, arc_y, 'Color',[0.7 0.1 0.4], 'LineWidth', 1.2);
    text(0.4 * r * cos(phi_rad/2), 0.4 * r * sin(phi_rad/2), '\phi','Color',[0.7 0.1 0.9], 'FontSize', 14, 'FontWeight', 'bold');

    % Title
    %title(sprintf('Lattice cell: a = %.2f, b = %.2f, phi = %.0f°', a, b, phi), 'FontSize', 12);
end
