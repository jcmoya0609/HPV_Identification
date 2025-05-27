%% Establishing variables and calculating MW and RSS
% EXAMPLE GRAIN ORIENTATION AND RS2C CHANGE THIS
v = [-0.406052; -0.053809; 0.129460]; %Rodrigues Vector orientation sample to crystal
RS2C = Rodrot(norm(v), unit(v)); %Ashley Bucsek's function to generate rotation matrix
cd /Users/janicemoya/Library/CloudStorage/OneDrive-Personal/Research/Matlab/Example_Data/
CTM = readmatrix("Shield_CTM_CuAlNi_Results.csv");

HPVNum = CTM(:,1); %EXAMPLE b AND m VALUES CHANGE THIS
b = CTM(:,4:6);
m = CTM(:,7:9); 
sigma=[0 0 0; 0 0 0; 0 0 10/0.0005^2*10^(-6)]; %CHANGE THIS 

[RSS_HPV_list, RSS_values, MW_HPV_list, MW_values] = HPVRankings(RS2C,m,b, sigma);
%% GrainMap Isosurface Render
% Example for GrainMapRender 
load cmapS1.mat
load grainsS1.mat  %EXAMPLE GRAIN MAP AND PARAMETERS
grainID=3;                  
parameters.patchmode = 'all';
parameters.cmap      = cmap; 
parameters.ambient   = 0.7;
parameters.diffuse   = 0.3;
parameters.specular  = 0.3;
parameters.lightaz   = -25;
parameters.lightel   = 0;
parameters.pxsize=0.0016;       %CHANGE THIS

% Call GrainMapRender to generate the isosurf data
[isosurf, parameters] = GrainMapRender(grainID, 0.45,cmap, parameters, grainsclean);

% ax = gca; % Get current axes
%     ax.XTick = []; % Remove tick marks for X-axis
%     ax.YTick = []; % Remove tick marks for Y-axis
%     ax.ZTick = []; % Remove tick marks for Z-axis


%% Plot HPV Plane
% Call HPVPlot function using the isosurf data

abc = [m(MW_HPV_list(1),:)]; % HPV plane coefficients  CHANGE number to whatever HPV you would like to test in the list of rankings by MW
z_offset2 = 0; % EXAMPLE Z_OFFSET, with this you can plot multiple parallel planes

region = struct('x_min', 0, 'x_max', 0.9, 'y_min', 0, 'y_max', 0.9, 'z_min', 0, 'z_max', 0.9); %CHANGE THIS TO YOUR 3D ROI WHERE YOUR GRAIN EXISTS
% Define plane color as a string (e.g., 'red', 'green')
planeColor = 'black'; % Color of the plane
% Call HPVPlot with the necessary inputs
if exist('h_patch', 'var') && ~isempty(h_patch)
    delete(h_patch(1)); %delete previously plotted plane if there exists one
else
    disp('h_patch does not exist or is empty.');
end
h_patch=HPVPlot(RS2C, abc, z_offset2, isosurf, grainID, parameters, region, planeColor);

