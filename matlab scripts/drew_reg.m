function drew_reg(filename, plotFlag)
    % drew_reg(filename, plotFlag)
    % plotFlag = 'scat' | 'sur'
    % for 'scatter' | 'surface'
    % Example:
    %   drew_reg('11_reg.dat','scat')
    %   drew_reg('11_reg.dat','con')
    %   drew_reg('11_reg.dat','sur')

    % Load points
    points = load(filename);

    % Separate into x, y, and value columns
    x = points(:,1);
    y = points(:,2);
    val = points(:,3);

    % Create grid for contour/surface if needed
    if ~strcmp(plotFlag,'scat')
        xlin = linspace(min(x), max(x), 100);
        ylin = linspace(min(y), max(y), 100);
        [X, Y] = meshgrid(xlin, ylin);
        Z = griddata(x, y, val, X, Y, 'cubic');
    end

    % Switch between plot types
    figure;
    switch lower(plotFlag)
        case 'scat'
            scatter(x, y, 50, val, 'filled');
            xlabel('X'); ylabel('Y');
            title(['Scatter plot from ', filename]);
            colorbar;

        case 'sur'
            surf(X, Y, Z);
            view(2);
            shading interp;
            colormap jet;
            xlabel('X'); ylabel('Y'); zlabel('Value');
            title(['Surface plot of ', filename]);
            
            colorbar;

        otherwise
            error('Unknown plotFlag. Use ''scatter'', ''contour'', or ''surface''.');
    end
end
