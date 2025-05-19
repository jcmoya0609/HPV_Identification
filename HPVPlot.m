function h_patch=HPVPlot(RS2C, abc, z_offset2, isosurf, ii, parameters, region, planeColor)
    % Function to plot the intersecting plane on the grain isosurface.
    % 
    % Inputs:
    % RS2C        - 3x3 matrix Sample to crystal grain averaged orientation
    % abc         - Coefficients of the HPV plane (1x3).
    % z_offset2   - Z-offset value to adjust the position of the plane.
    % isosurf     - Structure containing isosurface patch data.
    % ii          - The grain id index of the isosurface to plot.
    % parameters  - Parameters structure (contains pxsize, etc.).
    % region      - Structure containing x_min, x_max, y_min, y_max, z_min, z_max.
    % planeColor  - Color of the plane as a string (e.g., 'blue', 'red').

    % Perform the rotation matrix computation
    % R_S2C = Rodrot(norm(v), unit(v)); 
    
    a = abc(1); b = abc(2); c = abc(3);
    d = 0;  % Constant term
    
    % Define the range and resolution for the grid
    x_min = region.x_min; x_max = region.x_max;
    y_min = region.y_min; y_max = region.y_max;
    resolution = parameters.pxsize;  % Resolution of the grid

    iso_vertices1 = isosurf.patch{1};
    iso_vertices = iso_vertices1.vertices;

    % Extract x, y, z from iso_vertices
    x_new = iso_vertices(:, 1);
    y_new = iso_vertices(:, 2);
    z_new = iso_vertices(:, 3);

    % Calculate z values of the plane for the given (x, y)
    z_plane = (d - a * x_new - b * y_new) / c + z_offset2;

    % Combine plane coordinates and rotate them to be in sample frame
    crystal_coords = [x_new'; y_new'; z_plane'];  % 3 x N    
    sample_coords = transpose(RS2C) * crystal_coords;  % 3 x N in sample frame
    coord = sample_coords'; % Storing coordinates as N x 3

    % Check if each vertex is near the plane
    distances = abs(z_new - coord(:, 3));  % Difference in z-coordinates
    near_plane = distances < resolution;  % Logical array of points near the plane

    % Store visible coordinates
    visible_coords = iso_vertices(near_plane, :);  % Filter only near-plane vertices

    % Filter the visible coordinates based on the region limits
    filtered_coords = [];
    if ~isempty(visible_coords)
        filtered_coords = visible_coords(visible_coords(:,1) >= x_min & visible_coords(:,1) <= x_max & ...
                                         visible_coords(:,2) >= y_min & visible_coords(:,2) <= y_max & ...
                                         visible_coords(:,3) >= region.z_min & visible_coords(:,3) <= region.z_max, :);
    end

    % Perform Delaunay triangulation on the filtered coordinates
    if ~isempty(filtered_coords)
        tri_filtered = delaunay(filtered_coords(:,1), filtered_coords(:,2));

        % Plot the filtered intersecting patch with the specified plane color
        h_patch(1) = patch('Vertices', filtered_coords, 'Faces', tri_filtered, 'FaceColor', planeColor, ...
            'FaceAlpha', 0.75, 'EdgeColor', 'none');
        disp('done');
    else
        h_patch=[];
    end

    disp('completed');
end
