function cmap =CubeColormap(colormapName, numFaces)
    % APPLYCOLORMAPTOCUBE Generates a cube with colored faces using a colormap
    %
    %   cmap = applyColormapToCube(colormapName, numFaces)
    %
    % Inputs:
    %   colormapName - Name of the colormap (e.g., 'parula', 'hot', 'cool')
    %   numFaces - Number of faces in the cube (default: 6 for a regular cube)
    %
    % Output:
    %   cmap - The generated colormap used for the cube faces
    
    if nargin < 2
        numFaces = 6; % Default for a standard cube
    end
    
    % Define the cube vertices
    vertices = [
        0 0 0; % Vertex 1
        1 0 0; % Vertex 2
        1 1 0; % Vertex 3
        0 1 0; % Vertex 4
        0 0 1; % Vertex 5
        1 0 1; % Vertex 6
        1 1 1; % Vertex 7
        0 1 1; % Vertex 8
    ];

    % Predefined cube faces
    faces = [
        1 2 6 5; % Front face
        2 3 7 6; % Right face
        3 4 8 7; % Back face
        4 1 5 8; % Left face
        1 2 3 4; % Bottom face
        5 6 7 8; % Top face
    ];

    % Generate the colormap
    cmap = feval(colormapName, numFaces); % Use specified colormap
    
    % Create a figure for the cube
    figure;
    hold on;
    
    % Plot each face with its corresponding colormap color
    for i = 1:size(faces, 1)
        patch('Vertices', vertices, 'Faces', faces(i, :), ...
              'FaceColor', cmap(i, :), 'EdgeColor', 'k', 'LineWidth', 1.5);
    end
    
    % Configure plot appearance
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    view(3); % Set 3D view
    axis equal;
    grid on;
    hold off;
    
    % Return the colormap
    disp('Colormap generated and saved for reuse.');
end
