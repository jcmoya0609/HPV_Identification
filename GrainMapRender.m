

function [isosurf, parameters] = GrainMapRender(grainID,face_alpha, cmap, parameters, grainsclean)
    % grainNum:     - 1 x n array with Total number of grains in sample you need visualized eg 1 or [1:20]
    % face_alpha    - 1 x 1 double value between 0 and 1 representing the grain map transparency
    % cmap:         - n x 3 array with n being the number of grains in map. Assigns colormap to isosurface
    % parameters:   - 1x1 structure containing information about the lighting, pixel size and voxel size. Example: 
        % parameters.patchmode = 'all';
        % parameters.grainids= 1:grainNum;
        % parameters.cmap      = cmap; 
        % parameters.ambient   = 0.7;
        % parameters.diffuse   = 0.3;
        % parameters.specular  = 0.3;
        % parameters.lightaz   = -25;
        % parameters.lightel   = 0;
        % parameters.pxsize=0.0016;
    % grainsclean:   - m x n x p 3D array that contains the 3D
                    % Reconstruction of the grain map where each voxel is 
                    % assigned to a grainID

    %Establishing parameters structure
    parameters.grainids= grainID;
    parameters.cmap      = cmap;
    grainids=parameters.grainids; % If you would like to exclude certain grains you can change this

    % Initialize the isosurface structure
    isosurf = struct();
    isosurf.patch = cell(1, length(parameters.grainids));

    ll = figure;
    set(ll, 'Colormap', parameters.cmap);

    for ii = 1:length(parameters.grainids)
        current_grain_id = parameters.grainids(ii);

        if ismember(current_grain_id, grainids) 
            % face_alpha = 0.45;
            vol_bin    = (grainsclean == current_grain_id);
            vol_smooth = smooth3(vol_bin, 'box', 3);

            xx = (1:size(vol_smooth, 2)) * parameters.pxsize;
            yy = (1:size(vol_smooth, 1)) * parameters.pxsize;
            zz = (1:size(vol_smooth, 3)) * parameters.pxsize;

            % Store isosurface data
            isosurf.patch{ii} = isosurface(xx, yy, zz, vol_smooth, 0.5);
        
            % Render the patch
            h_patch(ii) = patch(isosurf.patch{ii}, 'EdgeColor', 'none', ...
                'FaceColor', parameters.cmap(current_grain_id, :), 'FaceAlpha', face_alpha);
        else
            face_alpha = 0; % If grain is excluded
        end
    end

    h_patch = findobj(gca, 'Type', 'patch');
    set(h_patch, 'AmbientStrength', parameters.ambient)
    set(h_patch, 'DiffuseStrength', parameters.diffuse)
    set(h_patch, 'SpecularStrength', parameters.specular)

    h_light = camlight(parameters.lightaz, parameters.lightel);

    % Customize axes appearance
     axis on
    view(3)
    axis equal

    % ax = gca; % Get current axes
    % ax.XTick = []; % Remove tick marks for X-axis
    % ax.YTick = []; % Remove tick marks for Y-axis
    % ax.ZTick = []; % Remove tick marks for Z-axis
    xlabel('X', 'FontSize', 18);
    ylabel('Y', 'FontSize', 18);
    zlabel('Z', 'FontSize', 18);
    set(gcf, 'Color', 'w'); % Change figure background to white
   
    

    disp('GrainMapRender completed');
end


           
   



