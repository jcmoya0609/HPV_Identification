% Written by Adam Creuziger (NIST)
% Modified by Ashley Bucsek and Celeste Perez (UM)
% Plots 2D and 3D EBSD and optical microscopy images of the CuAlNi sample from the ESRF ID06 June '22 beamtime
% Requires mtex version 6 or higher for 3 dimensional plotting of EBSD data (I used mtex v6.0beta2)
clc
% close all
clear

%% Questions

% 1. Diffty is in the middle - or does the scan start at diffty?


%% CHECK LOCAL PATH NAMES
% setMTEXpref('voronoiMethod','jcvoronoi');
mtexPath= '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/mtex-6.2.beta.3';  % Path to mtex folder
addpath(mtexPath);  startup_mtex
addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data');

%% Load previously analyzed EBSD data
load('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/PlottingIn3DScripts/S220607-AAC-009-rev1/S220607-AAC-009-rev1-2umStep-ws.mat')
addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/PlottingIn3DScripts/Mosa_ZSums_G7')
addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/PlottingIn3DScripts/Mosa_ZSums_G4')


%% USER INPUTS
% General plot settings for line and text sizes. Adjust to suit
ID_font_size = 24;
Grain_linewidth = 2;
Scan_linewidth = 4;

% crystal symmetry
CS = {'notIndexed', crystalSymmetry('m-3m', [5.8 5.8 5.8], 'mineral', 'color',...
    [0.53 0.81 0.98]),crystalSymmetry('mmm', [4.4 5.3 4.2], 'mineral', ...
    'CuAlNi-gammaprime', 'color', [0.56 0.74 0.56])};

% Set the X axis (i.e., the loading direction) for coloring
ipfKey = ipfHSVKey(ebsd_beta);
ipfKey.inversePoleFigureDirection = vector3d.Z;

% zsum layer stuff
fig_layers = 7;  % 6 = SMALL ROI OPTICAL MICROGRAPH, 7 = SMALL ROI EBSD, 8 = BIG OPTICAL MICROGRAPH, 9 = BIG EBSD

thresh = 130;  % Threshold choice for zsums (if want an additional threshold)
medFilterNeighb = [2 2];  % Median filter neighborhood size choice
% XSTART = 4905;  % THIS PARAMETER IS USED TO DECLARE THE FIRST LAYER'S VERTICAL POSITION -- SHIFTED BASED ON EYE!!!  --> There's some additional tweaking of the postions later too...see line 205 and lines 244-247
% LAYER_SPACING = 70;  % THIS IS THE LAYER SPACING IN MICRONS
xbins = 1 : 2560;  ybins = 1 : 2160*3;  % x and y - number of pixels  --> This is currently set up for 2X !!!
xbins = xbins * 208*1e-3;  ybins = ybins * 208*1e-3;  % x and y - convert from pixels to micrometers  --> This is currently set up for 2X !!!


%% Read in optical image data recorded at ESRF
I = imread("Creuziger-S220607-AAC-009-marks-assembled.png");

um_per_pixel=0.6321; % 1582 pix/mm
binning=1;

image_scale=binning*um_per_pixel;
I_scaled=imresize(I,image_scale);

% Set the scan width to draw representative lines
% From discussion with Celeste Perez
% 2X is 208nm per pixel
% 10X is 42nm per pixel
% 2560 pixel wide image

%%% -------ORIGINAL
% % 10X scan is 107 um wide
% scan_line_10X = 107/1000;%mm
% % 2X scan us 532 um wide
% scan_line_2X = 532/1000;%mm
% scan_line_choice = scan_line_10X;


%%% -------TEST
scan_line_10X = 107;   % µm
scan_line_2X  = 532;   % µm
scan_line_choice = scan_line_2X;


 %%
% transform orientations into a list of colors
color = ipfKey.orientation2color(ebsd_beta.orientations);

% ROI for close-up on grains
ROI_small=[3000 200 2000 1000];  % Create a smaller ROIROI_small=[3000 200 2000 1000];  % Create a smaller ROI
% ROI_small_g7=[4150 300 1000 800];  % Create a smaller ROI
ROI_smaller=[3300 500 500 300];  % Create a smaller ROIROI_small=[3000 200 2000 1000];  % Create a smaller ROI


ebsd_small = ebsd_beta(inpolygon(ebsd_beta,ROI_smaller));
gb_angle=5*degree;
% gb_angle=1*degree;

% Redoing the calc grains takes a while. Also, the grain IDs are not preserved, so can't use the numbering
[grains_beta_small,ebsd_small.grainId,ebsd_small.mis2mean] = calcGrains(ebsd_small,'angle',gb_angle,'boundary','tight');
color_small = ipfKey.orientation2color(ebsd_small.orientations);


%% PLOT SMALL ROI OPTICAL MICROGRAPH (3D) - VERTICAL - GRAYSCALE
% Per conversation at https://github.com/mtex-toolbox/mtex/issues/2014 and https://en.wikipedia.org/wiki/Isometric_projection

% '' Original - Don't change ''
ebsd_small.plottingConvention.north= vector3d.X;
ebsd_small.plottingConvention.outOfScreen=vector3d(1,-1,sqrt(2));
ebsd_small.plottingConvention.north= vector3d.X;
ebsd_small.plottingConvention.outOfScreen=vector3d(1,-1,sqrt(2));

% % '' Trying a different one ''
% ebsd_small.plottingConvention.north= vector3d.X;
% ebsd_small.plottingConvention.outOfScreen=vector3d(1,-0.5,2.3);
% ebsd_small.plottingConvention.north= vector3d.X;
% ebsd_small.plottingConvention.outOfScreen=vector3d(1,-0.5,2.3);
ebsd_small.plottingConvention
    % ans = plottingConvention
    % outOfScreen: (-1,-1,-1)
    % north : (0,1,-1)
    % east : (-2,1,1)

%% Fig 6 
close all
figure(6);
% figure;
plot(ebsd_small,ebsd_small.bc,ebsd_small.plottingConvention)%,'coordinates','on','labels','on')
colormap gray  % make the image grayscale
alpha(.2)  % Can change!
hold on
plot(grains_beta_small.boundary,'linewidth',Grain_linewidth)

z_ROI_front=0;
z_ROI_back = -138;

% %-------- original cropping 
% fill3( [5000 5000 5000 5000],...
%     [200 1200 1200 200],...
%     [z_ROI_front z_ROI_front z_ROI_back z_ROI_back],...
%     [0.3 0.3 0.3],'facealpha',0.4,'edgealpha',1)  % AB
% 
% fill3( [5000 5000 3000 3000],...
%     [200 200 200 200],...
%     [z_ROI_front z_ROI_back z_ROI_back z_ROI_front],...
%     [0.3 0.3 0.3],'facealpha',0.4,'edgealpha',1)  % AB

% 
%-------- smaller cropping G4
fill3( [3800 3800 3800 3800],...
    [500 800 800 500],...
    [z_ROI_front z_ROI_front z_ROI_back z_ROI_back],...
    [0.3 0.3 0.3],'facealpha',0.4,'edgealpha',1)  % AB

fill3( [3800 3800 3300 3300],...
    [500 500 500 500],...
    [z_ROI_front z_ROI_back z_ROI_back z_ROI_front],...
    [0.3 0.3 0.3],'facealpha',0.4,'edgealpha',1)  % AB
% % set(gcf, 'units', 'pixels', 'position', [-1899,276,477,788]);

return
%% COLORED (5) PLOT SMALL ROI EBSD (3D) - VERTICAL
figure(5);
plot(ebsd_small,color_small,ebsd_small.plottingConvention)%,'coordinates','on','labels','on')
alpha(.2)  % Can change!
hold on

z_ROI_front=0;
z_ROI_back = -138;

plot(grains_beta_small.boundary,'linewidth',Grain_linewidth)

fill3( [5000 5000 5000 5000],...
    [200 1200 1200 200],...
    [z_ROI_front z_ROI_front z_ROI_back z_ROI_back],...
    [0.3 0.3 0.3],'facealpha',0.4,'edgealpha',1)  % AB

fill3( [5000 5000 3000 3000],...
    [200 200 200 200],...
    [z_ROI_front z_ROI_back z_ROI_back z_ROI_front],...
    [0.3 0.3 0.3],'facealpha',0.4,'edgealpha',1)  % AB
% return
% view([90 95])
%% PLOT BIG OPTICAL MICROGRAPH (3D) - VERTICAL
% figure(8);
% plot(ebsd_beta,ebsd_beta.bc,ebsd_small.plottingConvention)%,'coordinates','on','labels','on')
% colormap gray  % Make the image grayscale
% alpha(.3)  % Can change!
% hold on
% plot(grains_beta.boundary,'linewidth',Grain_linewidth)
% 
% rectangle('position',ROI_small,'edgecolor','w','linewidth',2)  % Draw ROI box if desired
% 
% fill3( [6.4986e+03 6.4986e+03 6.4986e+03 6.4986e+03],...
%     [0 1.4100e+03 1.4100e+03 0],...
%     [z_ROI_front z_ROI_front z_ROI_back z_ROI_back],...
%     [0.3 0.3 0.3],'facealpha',0.3,'edgealpha',1)  % AB
% 
% fill3( [6.4986e+03 6.4986e+03 0 0],...
%     [0 0 0 0],...
%     [z_ROI_front z_ROI_back z_ROI_back z_ROI_front],...
%     [0.3 0.3 0.3],'facealpha',0.3,'edgealpha',1)  % AB
% 
% text(grains_beta,grains_beta.id, 'FontSize', 24)

%% PLOT GRAIN 7 - BIG EBSD (3D) - VERTICAL
% figure(9);
% plot(ebsd_beta,color,ebsd_small.plottingConvention)%,'coordinates','on','labels','on')
% alpha(.5)  % Can change!
% hold on
% plot(grains_beta.boundary,'linewidth',Grain_linewidth)
% 
% rectangle('position',ROI_small,'edgecolor','w','linewidth',2)  % Draw ROI box if desired
% 
% fill3( [6.4986e+03 6.4986e+03 6.4986e+03 6.4986e+03],...
%     [0 1.4100e+03 1.4100e+03 0],...
%     [z_ROI_front z_ROI_front z_ROI_back z_ROI_back],...
%     [0.3 0.3 0.3],'facealpha',0.3,'edgealpha',1)  % AB
% 
% fill3( [6.4986e+03 6.4986e+03 0 0],...
%     [0 0 0 0],...
%     [z_ROI_front z_ROI_back z_ROI_back z_ROI_front],...
%     [0.3 0.3 0.3],'facealpha',0.3,'edgealpha',1)  % AB
% 
% text(grains_beta,grains_beta.id, 'FontSize', 24)


%% GRAIN 7 - VERTICAL


fig_layers = 6;

matfile_name_G7 = '9p0_BandW_mosalayer_2x_%02d.mat';
% matfile_name_G4 = 'G7_8p3_Brown_mosalayer_2x_%02d.mat';
XSTART = 4920;  % THIS PARAMETER IS USED TO DECLARE THE FIRST LAYER'S VERTICAL POSITION -- SHIFTED BASED ON EYE!!!  --> There's some additional tweaking of the postions later too...see line 205 and lines 244-247
LAYER_SPACING = 50;  % THIS IS THE LAYER SPACING IN MICRONS
range = [1 2 9];
for i = 1:length(range)

layernum = range(i);
fprintf('Processing Layer: %.0f ...\n', layernum);    
filename = sprintf(matfile_name_G7, layernum);
load(filename);
zSumm = sprintf('zSum_%02d', layernum);
zSumm = BandW;zSum_01_unedited = BandW;
% zSumm =rot90(zSumm,4);
scalingFactor = mean(mean(zSum_01_unedited(~isnan(zSum_01_unedited)))) / mean(mean(zSumm(~isnan(zSumm))));
zSumm = zSumm * scalingFactor;  % Increase intensity to match layer 1
zSumm(zSumm<thresh) = nan;  % Threshold image
zSumm = medfilt2(zSumm, medFilterNeighb);  % Apply median filter
figure('Visible', 'off'); imagesc(xbins,ybins,log(zSumm)); colormap bone; colorbar; axis square; axis equal;
set(gca,'DataAspectRatio',[1 1 1]); %ylim([2800 4000])
diffty = 0.05; difftz =-LAYER_SPACING/100;  % RANDOMLY SELECTED
[x_ebsd, y_ebsd]=ESRF2EBSD(diffty,difftz);
zPos = x_ebsd;
clearvars h h2 X Y Z C xx yy zz 
[xx,yy,zz] = meshgrid(xbins,ybins,zPos);

figure('Visible', 'off'); h = surf(xx,yy,zz,log(zSumm),'edgecolor','none');
axis square;  axis equal; colormap bone; box on;  grid off;
set(gca,'DataAspectRatio',[1 1 1])
axis on; xlabel('x'); ylabel('y'); zlabel('z'); %view([10 15])

rotate(h, [0 1 0], 90.7071, [0 10 0]); %chi
rotate(h,[1 0 0],90,[0 10 0]);
rotate(h, [0 0 1], 0.3878, [0 00 0]);

axis vis3d equal
X = get(h,'xdata');  X = X - 150 ;  
X(:) = XSTART - LAYER_SPACING * (layernum);

Y = get(h,'ydata'); miny = min(Y);
Y = get(h,'ydata');  Y = Y + (y_ebsd - min(Y,[],'all') - scan_line_choice/2) ;
Z = get(h,'zdata');   %Z = flipud(Z); 
Z = Z-300;
% Z = Z - 65 - 133;
C = get(h,'cdata');
figure(fig_layers);lims = clim;
slopey = (lims(1)-lims(2))/((min(min(C)))-max(max(C)));
yinty = lims(1) - slopey*(min(min(C)));
C = C*slopey+yinty;
% bone_map = bone(256);
% bone_map = bone_map .^ 0.6;  % Make the colormap darker
% colormap(bone_map);
colormap(bone)
clim ([50 350])
figure(fig_layers); hold on;  h2 = surf(X,Y,Z,C,'edgecolor','none'); 
zlim([-138.5 0.5])
fprintf('Finished Processing Layer: %.0f ...\n', layernum);    
set(gca,'ytick',[]); set(gca,'ztick',[]);

end

% xlim padded
% end
% return
%% PLOT Grain 4 - VERTICAL loop
fig_layers = 6;
addpath /Volumes/MyPassport/ESRF2022_ID06/Celeste_mat_files/G4/
thresh = 70;
medFilterNeighb = [2 2];
xbins = (1:2560) * 42e-3;
ybins = (1:(3*2160)) * 42e-3;
XSTART = 3760-70;
LAYER_SPACING = 70;
y_extra = 60
Diffry = 5.97025;  %#ok<NASGU>  % kept for clarity
Chi    = 0.6051;   %#ok<NASGU>
flipbeam = 0;

%%Reference mean from layer 00 (UNEDITED)
load 8p0_2_mosalayer_10x_00.mat
zSum_01_unedited = inTot;                             % ref (unedited)
refMean = mean(zSum_01_unedited(~isnan(zSum_01_unedited)), 'all');

% %%Loop over layers 00..04 (flip on 00–02 only)
% files = { '8p0_2_mosalayer_10x_00.mat', ...
%           '8p0_2_mosalayer_10x_01.mat', ...
%           '8p0_2_mosalayer_10x_02.mat', ...
%           '8p0_2_mosalayer_10x_03.mat', ...
%           '8p0_2_mosalayer_10x_04.mat' };
% % flipLayer = [true true true true true];  %
% flipLayer = [false false false false false];  %

% ---- file list via formatted name ----
matfile_name_G7 = '8p0_2_mosalayer_10x_%02d.mat';
range =1:4;                                % layers 00..04
flipLayer = false(1, numel(range));         % set true/false per layer as needed

for k = 1:numel(range)
    layernum = range(k);
    fprintf('Processing Layer: %02d ...\n', layernum);
    filename = sprintf(matfile_name_G7, layernum);

    % SAFE LOAD
    S = load(filename, 'inTot');            % load just inTot from file
    assert(isfield(S,'inTot'), 'File %s missing variable inTot.', filename);
    z = S.inTot;

    if flipLayer(k)
        z = flip(z, 1);
    end
    % Scale to layer00 mean
    curMean = mean(z(~isnan(z)), 'all');
    z = z * (refMean / curMean);

    % Threshold + median filter
    z(z < thresh) = NaN;
    z = medfilt2(z, medFilterNeighb);

    % ---- EBSD mapping & gri90 ----
    diffty = 0.05;  difftz = -0.7;
    [x_ebsd, y_ebsd] = ESRF2EBSD(diffty, difftz);
    zPos = x_ebsd;
    [xx,yy,zz] = meshgrid(xbins, ybins, zPos);

    % ---- Build rotated surface on a temp (invisible) fig, then extract ----
    ftmp = figure('Visible','off');
    h = surf(xx,yy,zz, log(z), 'edgecolor','none');

    rotate(h,[0 1 0], 96.1164,        [0 0 0]);   % chi
        % rotate(h,[0 1 0], 120.1164,        [0 0 0]);   % chi

    rotate(h,[1 0 0], 90,             [0 0 0]);
    rotate(h,[0 0 1], flipbeam+0.765, [0 0 0]);   % diffry
axis vis3d equal
    X = get(h,'xdata');  X = X - 150;  X(:) = XSTART - LAYER_SPACING*(k);
    Y = get(h,'ydata');  Y = Y + (y_ebsd - min(Y,[],'all') - scan_line_10X/2) - y_extra;
    Z = get(h,'zdata');  Z = Z - 65 - 133;  Z = flipud(Z);
    C = get(h,'cdata');
    close(ftmp);

    % ---- Color remap to current figure limits, then plot there ----
    figure(fig_layers); lims = clim;           % use existing clim on fig 6
    slopey = (lims(1)-lims(2)) / (min(C,[],'all') - max(C,[],'all'));
    yinty  = lims(1) - slopey * min(C,[],'all');
    C = C * slopey + yinty;
colormap(bone); clim ([0 2600])
    figure(fig_layers); hold on
    surf(X, Y, Z, C, 'edgecolor','none');
    axis vis3d equal




end

% Optional: tidy ticks (as in your last lines)
set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'ztick',[]);

fprintf('Finished Processing Layer');    
view([-90 74]) % VIEWWWW



return



%% Grain 7 corrected

% User Inputs
fig_layers       = 6;
matfile_name_G7  = '9p0_BandW_mosalayer_2x_%02d.mat';
Diffry           = 0.37215;  
Chi              = 0.7071; 
diffty           = 0.05; 
difftz           = -0.5; 
XSTART           = 4920;   % FIRST LAYER'S VERTICAL POSITION 
LAYER_SPACING    = 50;     % LAYER SPACING IN MICRONS
thresh = 70;
range            = [1 2 9];
scan_line_choice = scan_line_2X;

view([-90 74]) % VIEW

    

% Loop Starts:
for i = 1:length(range)
    layernum = range(i);
    fprintf('G7 - Processing Layer: %.0f ...\n', layernum);
    filename = sprintf(matfile_name_G7, layernum);
    load(filename);

    zSumm               = sprintf('zSum_%02d', layernum);
    zSumm               = BandW; % .mat variable name
    zSum_01_unedited    = BandW; % .mat variable name
   

    % Scaling Factor: Normalize brightness to match reference
    scalingFactor = mean(mean(zSum_01_unedited(~isnan(zSum_01_unedited)))) / mean(mean(zSumm(~isnan(zSumm)))); 
    zSumm = zSumm * scalingFactor; 


    % Threshold image: 
    zSumm(zSumm<thresh) = nan; 


    % Apply median filter:
    zSumm = medfilt2(zSumm, medFilterNeighb);  
    set(gca,'DataAspectRatio',[1 1 1]); 
   

    % Computes where this layer should sit in EBSD coordinates (ESRF2EBSD)
    [x_ebsd, y_ebsd]=ESRF2EBSD(diffty,difftz); 
    zPos = x_ebsd; 
    clearvars h h2 X Y Z C xx yy zz 


    % Uses a meshgrid to map xbins, ybins into 3D space at zPos
    [xx,yy,zz] = meshgrid(xbins,ybins,zPos); 


    % Creates a temporary 3D surface from the intensity image
    figure('Visible', 'off'); 
    h = surf(xx,yy,zz,log(zSumm),'edgecolor','k'); axis square; axis equal; 
    colormap bone; box on; grid off; set(gca,'DataAspectRatio',[1 1 1]); axis on; 


    %Apply motor rotations
    Chi7 = 90 + Chi;
    % Diffry7 = 90 + Diffry;
    % rotate(h, [0 1 0], Chi7 , [0 0 0]); % Chi
    % rotate(h,[1 0 0],90,[0 0 0]);
    % rotate(h, [0 0 1], Diffry7, [0 0 0]);  % Diffry
    %Apply motor rotations
    rotate(h, [0 1 0], Chi7 , [0 0 0]); % Chi
    rotate(h,[1 0 0],90,[0 0 0]);
    rotate(h, [0 0 1], Diffry, [0 00 0]); % Diffry

    axis vis3d equal % Makes sure no distorsion
    

    % Sets the X-position based on layer number (stacking layers).
    X = get(h,'xdata'); 
    X = X - 150; 
    X(:) = XSTART - LAYER_SPACING * (layernum - 1); 


    % Adjusts Y for EBSD alignment.
    Y = get(h,'ydata'); 
    Y = Y + (y_ebsd-min(min(Y))-scan_line_choice/2); 


    % Shifts Z to align with the global frame.
    Z = get(h,'zdata'); Z = Z-590-133; %Z = flipud(Z); 


    % Rescales intensity values (C) so all layers share the same colormap scale.
    C = get(h,'cdata'); 
    figure(fig_layers); lims = clim; 
    slopey = (lims(1)-lims(2))/((min(min(C)))-max(max(C))); 
    yinty = lims(1) - slopey*(min(min(C))); 
    C = C*slopey+yinty; 


    % Final Plot into Fig_layers
    figure(fig_layers);
    hold on; h2 = surf(X,Y,Z,C,'edgecolor','none'); 
    zlim([-138.5 0.5]) 
    fprintf('G7 - Finished processing Layer: %.0f ...\n', layernum);
end
%% Grain 4 corrected
% fig_layers = 6;
% addpath /Volumes/MyPassport/ESRF2022_ID06/Celeste_mat_files/G4
% addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/ThreeD_matfigures/Grain_4/Mosa_ZSums_G4/')
addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/ThreeD_matfigures/Grain_4/Mosas/8p0_Mosas/')

xbins = (1:2560) * 42e-3;
ybins = (1:(3*2160)) * 42e-3;
% XSTART = 3760;
% LAYER_SPACING = 70;
% y_extra = 60
% Diffry = 6.1164;  %#ok<NASGU>  % kept for clarity
% Chi    = 0.765;   %#ok<NASGU>
% flipbeam = 0;



% User Inputs
fig_layers       = 6;
% matfile_name_G4  = '8p0_Brown_mosalayer_10x_%02d.mat';
matfile_name_G4  = '8p0__ZSUM_Mosas_10x_%02d.mat';
Diffry           = 5.97025;  
Chi              = 0.6051; 
diffty           = 0.05; 
difftz           = -0.43; 
XSTART           = 3700;   % FIRST LAYER'S VERTICAL POSITION 
LAYER_SPACING    = 70;     % LAYER SPACING IN MICRONS
thresh           = 70;
range            = [1 2 3];
scan_line_choice = scan_line_10X;

view([-90 74]) % VIEW
% xbins = 1 : 2560;  ybins = 1 : 3*2160;  % x and y - number of pixels
% xbins = xbins * 42*1e-3;  ybins = ybins * 42*1e-3;  % x and y - convert from pixels to micrometers

% Loop Starts:
for i = 1:length(range)
    layernum = range(i);
    fprintf('G4 Processing Layer: %.0f ...\n', layernum);
    filename = sprintf(matfile_name_G4, layernum);
    load(filename);
    % clearvars zSumm zSum_01_unedited
    zSumm               = sprintf('zSum_%02d', layernum);
    zSumm               = Zsumm; % .mat variable name
    zSum_01_unedited    = Zsumm; % .mat variable name
   

    % Scaling Factor: Normalize brightness to match reference
    scalingFactor = mean(mean(zSum_01_unedited(~isnan(zSum_01_unedited)))) / mean(mean(zSumm(~isnan(zSumm)))); 
    zSumm = zSumm * scalingFactor; 


    % Threshold image: 
    zSumm(zSumm<thresh) = nan; 


    % Apply median filter:
    zSumm = medfilt2(zSumm, medFilterNeighb);  
    set(gca,'DataAspectRatio',[1 1 1]); 
   
    % figure; imagesc(xbins,ybins,log(zSumm)); colormap bone; colorbar;
    % axis square; % axis equal; set(gca,'DataAspectRatio',[1 1 1]);
    % ylim([2800 4000])

    % Computes where this layer should sit in EBSD coordinates (ESRF2EBSD)
    [x_ebsd, y_ebsd]=ESRF2EBSD(diffty,difftz); 
    zPos = x_ebsd; 
    clearvars h h2 X Y Z C xx yy zz 


    % Uses a meshgrid to map xbins, ybins into 3D space at zPos
    [xx,yy,zz] = meshgrid(xbins,ybins,zPos); 


    % Creates a temporary 3D surface from the intensity image
    figure('Visible', 'off'); 
    % figure
    h = surf(xx,yy,zz,log(zSumm),'edgecolor','k'); axis square; axis equal; 
    colormap bone; box on; grid off; set(gca,'DataAspectRatio',[1 1 1]); axis on; 
    
    
    %Apply motor rotations
    Chi0 = 90 + Chi;
    rotate(h, [0 1 0], Chi0 , [0 0 0]); % Chi
    rotate(h,[1 0 0],90,[0 0 0]);
    rotate(h, [0 0 1], Diffry, [0 0 0]);  % Diffry

    axis vis3d equal % Makes sure no distorsion
    

    % Sets the X-position based on layer number (stacking layers).
    X = get(h,'xdata'); 
    X = X - 150; 
    X(:) = XSTART - LAYER_SPACING * (layernum - 1); 


    % Adjusts Y for EBSD alignment.
    Y = get(h,'ydata'); 
    Y = Y + (y_ebsd-min(min(Y))-scan_line_choice/2); 


    % Shifts Z to align with the global frame.
    Z = get(h,'zdata'); Z = Z-65-140; Z = flipud(Z); 


    % Rescales intensity values (C) so all layers share the same colormap scale.
    C = get(h,'cdata'); 
    figure(fig_layers); lims = clim; 
    slopey = (lims(1)-lims(2))/((min(min(C)))-max(max(C))); 
    yinty = lims(1) - slopey*(min(min(C))); 
    C = C*slopey+yinty; 


    % Final Plot into Fig_layers
    figure(fig_layers);
    hold on; h2 = surf(X,Y,Z,C,'edgecolor','none'); 
    % zlim([-138.5 0.5]) 
    fprintf('G4 Finished processing Layer: %.0f ...\n', layernum);
end



%% Mixing codes 
mtexPath= '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/mtex-6.2.beta.3';  % Path to mtex folder
addpath(mtexPath);  startup_mtex
setMTEXpref('voronoiMethod','jcvoronoi');setMTEXpref('voronoiMethod','jcvoronoi');





%% Plots individual crystals with variants from fig 81 - works
idx=7; % Grain 4 with prior labelling
% idx=9; % Grain 5 with prior labelling
% idx=3; % Grain 7 with prior labelling

ori=grains.meanOrientation(idx);
%%Draw individual grains
figure(51)
plot(cSGrains(idx),'coordinates','on','faceAlpha',0.5);
axis on; 
view(0,90);

range = numVectors;
% range = [1 54 60 80];
% Loop through all vectors in sS
for n = 1:length(range)
    k = range(n)
    figure(100+k)
    % figure('Name','Variant #',k)

    % Plot the crystal shape
    plot(ori * cS, 'faceAlpha','facecolor',[0.1216,0.9804 , 0.3098], 0.2, 'LineWidth', 1);
    
    % Title for each plot
    title(['\textbf{' int2str(k) '}:' char(sS(k).n, 'latex')], ...
          'Interpreter', 'latex');
    axis on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    

    hold on;
    
    % Plot the rotated shape or vector with a color from the colormap
    plot(ori * cS, ori * sS(k), 'facecolor', cmap(k, :),0.8, ...
         'arrowLineWidth', 3, 'LineWidth', 1.5);
    
    % % Add an arrow (customizable appearance)
    % arrow3d(0.4 * xvector, 'faceColor', 'k', 'linewidth', 3);
    
    hold off;
    view(0,90);
end
