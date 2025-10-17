%% Built from of CuAlNi_AC_CP_s5_b1.m
% Commit 67f8be50fc22ec2a816eb0a5f8a1f586fef21f79
% Branch "EBSD-import"

% Plotting wireframe crystal orientation
% Plotting Austenite-Twinned Marteniste (m) vectors
% Plotting rotated m vectors
% Using slip & schmid functions in mtex to do avail work

warning('off')

%% CLEAR ALL PRIOR OUTPUT
clc; clear; close all;
mtexPath= '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/mtex-6.2.beta.3';  % Path to mtex folder
addpath(mtexPath);  startup_mtex


%% SET-UP PATH
% only needed initialy

setMTEXpref('voronoiMethod','jcvoronoi');
%mtexPath= '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/mtex-6.2.beta.3';  % Path to mtex folder
%addpath(mtexPath);  startup_mtex
%addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data');

%% CRYSTAL AND SPECIMEN AND SYMMETRIES 

% CRYSTAL SYMMETRY

CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [5.8 5.8 5.8], 'mineral',...
  'CuAlNi-beta', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('mmm', [4.4 5.3 4.2], 'mineral',...
  'CuAlNi-gammaprime', 'color', [0.56 0.74 0.56])};

% PLOTTING CONVENTION
setMTEXpref('xAxisDirection','east');
setMTEXpref('yAxisDirection','north');
%setMTEXpref('zAxisDirection','IntoPlane');



%% SPECIFY FILE NAMES 

%RELATIVE PATH TO FILES

pname = './Example_Data'; % Adam
%pname = '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data'; % Celeste
fname = [pname '/E220614-AAC-009_2umstep.ctf'];
% SAVE NAME AND PATH
plotname='CuAlNi_009';
savepath='CuAlNi_AC_009';
phase_name='CuAlNi-beta';

%% IMPORT THE DATA 

ebsd_org = EBSD.load(fname,CS,'interface','ctf',...
    'convertEuler2SpatialReferenceFrame','setting 2');

% Sample alignment determined by checking Oxford data and
% EBSD wireframe overlay
% Need to change the plot image, but not the angles
% Flip vertically, but orientations are ok
rot = rotation.byAxisAngle(xvector,180*degree);
%ebsd_org = rotate(ebsd_org,rot,'keepEuler');
ebsd_org = rotate(ebsd_org,rot);


% Need to change the angles, but keep the EBSD XY 
%rot = rotation.byAxisAngle(yvector,0*degree);
%ebsd_org = rotate(ebsd_org,rot,'keepXY');

% Rotate the EBSD data to align with sample edge
rotation_angle=2*degree
% define a rotation
rot = rotation.byAxisAngle(zvector,rotation_angle);

% rotate the EBSD data
ebsd_org = rotate(ebsd_org,rot);

%% Crop the data to reduced area

ROI = [600 -1550 6500 1300];
ROI_over = [350 -1600 7000 1400]; 

ebsd=ebsd_org(inpolygon(ebsd_org,ROI));

% change plotting convention, don't rotate the data...
%ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,180*degree))

%% FLAG TO SAVE FILES AND NAME 

%save_files=true;
save_files=false;

if not(isfolder(savepath))
    mkdir(savepath);
end

disp('####################################')
disp(strcat('Sample ', plotname))
disp('####################################')


%% From Janice's Code:

CTM = readmatrix(strcat(pname, "/Shield_CTM_CuAlNi_Results.csv"));
number_of_interfaces=96;
HPVNum = CTM(:,1);
b = CTM(:,4:6);
m = CTM(:,7:9); 
m_3dvec = vector3d(m');
b_3dvec = vector3d(b');
m_double=reshape(double(m_3dvec),[number_of_interfaces,3]); %csv
b_double=reshape(double(b_3dvec),[number_of_interfaces,3]); %csv 

%%=======Orientations from Shield 1995
% Mtex expects these as a Z and X pair
% Images and intercepts likely need to be rotated

A_x=orientation.byMiller([-0.380  0.925 0],[ 0.925  0.380 0 ],CS{2});
A1_T0=A_x;
figure(100);

h_JM = Miller({1,0,0},{1,1,0},{1,1,1}, CS{2});
plotPDF(A1_T0, h_JM, 'antipodal', 'MarkerSize',15,'marker','s',...
    'MarkerEdgeColor','r','MarkerFaceColor','r' );
hold on

%% ========================== Figures Start ===============================
%% PLOT BAND CONTRAST 

ebsd.how2plot.east = xvector;
ebsd.how2plot.north = yvector;
ebsd.how2plot
%% Fig (1) - Optical image 

fprintf('Processing figure (1) - Band Contrast\n');
figure(1);
%plot(ebsd,ebsd.bc,'micronbar','off');
plot(ebsd_org,ebsd_org.bc,'coordinates','on');
colormap gray; % this makes the image grayscale
mtexColorbar;
rectangle('position',ROI_over,'edgecolor','w','linewidth',2);
rectangle('position',ROI,'edgecolor','y','linewidth',2);
if save_files==true
    saveFigure(strcat(savepath, plotname, '-FullBandContrast.png'));
end
%% Fig (2) - EXPLICIT ORIENTATION COLORING
fprintf('Processing figure (2) - IPF Map\n');

% % DEFINA AN EXPLICIT MAP FOR CONVERING ORIENTATIONS INTO COLORS
oM = ipfHSVKey(CS{2});
% % CHANGE TO NEW DIRECTION:
% oM.inversePoleFigureDirection = zvector;

figure(2);
plot(oM);
if save_files==true
    saveFigure(strcat(savepath,plotname, '-IPF-Key.png'));
end

%% Fig (3) - POLE FIGURES
% Maybe do this after finding grain boundaries to limit the number of
% points.

% ORIENTATION AND MILLER INDICES
ebsd_phase = ebsd(phase_name)
h = Miller({1,0,0},{1,1,0},{1,1,1}, CS{2}); % 100 - face ; 110 edge ; 111 corner
ori = ebsd_phase.orientations

% POLE FIGURE

figure (3);
plotPDF(ori, h, 'antipodal', 'MarkerSize', 4);
title('Pole Figure - CuAlNi-beta');

% return
%% Plot UNIT CELL
% General Shape of crystal
figure(4);
cS = crystalShape.cube(ebsd_phase.CS)
plot(cS,'faceAlpha',0.2); hold on; 
%% Crop and Resample?
% could crop and resample to reduce the size of the data...

%% Fig (11, 12, 13, 14) - IPF MAPS

%%================ IPF X ================= %%
fprintf('Processing figure (11) - EBSDX\n');
oM.inversePoleFigureDirection = xvector;
color_unfiltered = oM.orientation2color(ebsd(phase_name).orientations);
figure(11);
plot(ebsd(phase_name),color_unfiltered,'figSize','large');
if save_files==true
    saveFigure(strcat(savepath,plotname, '-IPF-EBSDX-Raw.png'));
end

%%================ IPF Y ================= %%
fprintf('Processing figure (12) - EBSDY\n');
oM.inversePoleFigureDirection = yvector;
color_unfiltered = oM.orientation2color(ebsd(phase_name).orientations);
figure(12);
plot(ebsd(phase_name),color_unfiltered,'figSize','large');
if save_files==true
    saveFigure(strcat(savepath,plotname, '-IPF-EBSDY-Raw.png'));
end

%%================ IPF Z ================= %%
fprintf('Processing figure (13) - EBSDZ\n');
oM.inversePoleFigureDirection = zvector;
color_unfiltered = oM.orientation2color(ebsd(phase_name).orientations);
figure(13);
plot(ebsd(phase_name),color_unfiltered,'figSize','large');
if save_files==true
    saveFigure(strcat(savepath,plotname, '-IPF-EBSDZ-Raw.png'));
end
%% ================ IPF Maps Subplot ================= %%
fprintf('EBSD Subplot\n'); 

figure(14)
mtexFig = newMtexFigure('layout', [2, 2]);
phase_name_x = phase_name;  phase_name_y = phase_name;  phase_name_z = phase_name;
% IPF-X
oM_x = ipfHSVKey(ebsd(phase_name));  oM_x.inversePoleFigureDirection = xvector;
color_x = oM_x.orientation2color(ebsd(phase_name_x).orientations);
plot(ebsd(phase_name_x), color_x);  title('IPF-X');  nextAxis;
% IPF-Y
oM_y = ipfHSVKey(ebsd(phase_name));  oM_y.inversePoleFigureDirection = yvector;
color_y = oM_y.orientation2color(ebsd(phase_name_y).orientations);
plot(ebsd(phase_name_y), color_y);  title('IPF-Y');  nextAxis;
% IPF-Z
oM_z = ipfHSVKey(ebsd(phase_name));  oM_z.inversePoleFigureDirection = zvector;
color_z = oM_z.orientation2color(ebsd(phase_name_z).orientations);
plot(ebsd(phase_name_z), color_z);  title('IPF-Z');  nextAxis;
% IPF-Key
plot(oM);  title('IPF-Key');
% Add a figure-wide title
annotation('textbox', [0.6598 0.3847 0.0922 0.053],'String', '$m\bar{3}m$', ...
'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 16, ...
'FontWeight', 'bold','Interpreter', 'latex');
sgtitle('EBSD - IPF Maps','FontSize', 18,'FontWeight', 'bold');

if save_files==true
    saveFigure(strcat(savepath,plotname, '-IPF-Combined-Raw.png'));
end

%% SAVE ORIGINAL EBSD DATA BEFORE FILTERING
% unneeded, just adds storage

% ebsd_orig = ebsd; % Grain Reconstruction
% ebsd_unfiltered = ebsd; % Start from original
% ebsd_filled = ebsd;
% ebsd_roi = ebsd; 
% ebsd_color = ebsd;
% 
% % COLOR FROM ORIGINAL DATA UNFILTERED
% color_unfiltered = oM.orientation2color(ebsd_color(phase_name).orientations);

%% GRAIN RECONSTRUCTION

% INITIAL GRAIN RECONSTRUCTION
[grains, ebsd(phase_name).grainId] = calcGrains(ebsd(phase_name), 'angle', 3*degree);
% REMOVE SMALL GRAINS
ebsd(grains(grains.numPixel < 1000)) = [];
% RECOMPUTE GRAINS AFTER REMOVAL
[grains, ebsd(phase_name).grainId] = calcGrains(ebsd(phase_name), 'angle', 3*degree);
% SMOOTH THE GRAIN BOUNDARIES (OPTIONAL)
grains = smooth(grains, 5);

%% Could remap the Grain IDs

%% Di
% figure(21);
% Display average orientation of remaining grains in text
% disp(grains.meanOrientation);

%% Fig (22) -  EBSD + GRAIN BOUNDARIES (CP)
% NUMBERING IS DIFFERENT THAN WHEN ORIGINALLY PLOTTED
% likely a change from v5.x to v6.x
fprintf('Processing figure (22) -  EBSD + GRAIN BOUNDARIES (CP)\n');
% PLOT EBSD UNFILTERED WITH GB 
figure(22); clf;

color_x = oM_x.orientation2color(ebsd(phase_name).orientations);
plot(ebsd(phase_name), color_x, 'figSize', 'huge','coordinates','on');
hold on;
% PLOT ORIENTATION ON TOP
%plot(ebsd(phase_name), ebsd(phase_name).orientations,...
%    'figSize', 'large');
% PLOT GRAIN BOUNDARIES FROM FILTERED DATA
plot(grains.boundary, 'linewidth', 2);
text(grains,grains.id, 'FontSize',20)
hold off;

%%=============  Save Figure if Needed =============== %%
if save_files == true
    saveFigure(strcat(savepath, plotname, '-IPF-EBSDZ-3degGB-Raw.png'));
end


%% Fig (23) -  EBSD + GRAIN BOUNDARIES (CP) + Crystal Shape
fprintf('Processing figure (23) -  EBSD + GRAIN BOUNDARIES (CP) + Crystal Shape\n');

% SELECT IPF COLOR
% color_z = oM_z.orientation2color(ebsd(phase_name).orientations);
color_x = oM_x.orientation2color(ebsd(phase_name).orientations);

% PLOT EBSD UNFILTERED WITH GB 
figure(23);
cSGrains = grains.meanOrientation * cS * 0.7* sqrt(grains.area);
ebsd_color = color_x;
% grain_color = ebsd_color(grains(i).meanOrientation == ebsd(phase_name).orientations, :);
plot(ebsd(phase_name), ebsd_color, 'figSize', 'huge','coordinates','on');
hold on;
% PLOT ORIENTATION ON TOP
%plot(ebsd(phase_name), ebsd(phase_name).orientations,...
%    'figSize', 'large');
% PLOT GRAIN BOUNDARIES FROM FILTERED DATA
plot(grains.boundary, 'linewidth', 2);
text(grains,grains.id, 'FontSize',20);
% plot Crystal Shape

% % % CP - Adjust the centroids with an offset
% plot(grains(i).centroid + cSGrains(i), ...
%      'faceColor', grain_color, 'edgeColor', 'k', ...
%      'LineWidth', 2, 'faceAlpha', 0.7);

% % ORIGINAL ADAM
plot((grains.centroid+ cSGrains), 'faceColor', [1, 0.6, 0.6], ...
    'edgeColor', 'k', 'LineWidth', 2, 'faceAlpha', 0.7);
hold off;

%%=============  Save Figure if Needed =============== %%
if save_files == true
    saveFigure(strcat(savepath, plotname, '-IPF-EBSDZ-3degGB-Raw.png'));
end


%% Plot m vectors
figure(41)
scatter(m_3dvec,'grid','on','antipodal')

figure(42)
scatter(m_3dvec,'grid','on')

%% 51 - Choose a grain
idx=7; % Grain 4 with prior labelling
% idx=9; % Grain 5 with prior labelling
% idx=3; % Grain 7 with prior labelling

ori=grains.meanOrientation(idx);

%%Draw individual grains

figure(51)
plot(cSGrains(idx),'coordinates','on','faceAlpha',0.5);
axis on; hold on
xlabel X; ylabel Y; zlabel Z;
hold off
% Set the 3D view (azimuth, elevation)
% (0, 270) gives the X-Y plane with Y positive downward
% need to match what the EBSD plotting axes display
view(0,90);
% view(3)


%% Plot m vectors rotated

%WITH ROTATION 
figure(61)
scatter(grains.meanOrientation(idx)*m_3dvec,'grid','on','antipodal');


%% 80 - Define using slipSystem in mtex W/O rotations
fprintf('Processing figure (80) -  Define using slipSystem in mtex W/O rotations\n');

sS = slipSystem(b_3dvec, m_3dvec);

figure(80);
% For some reason, need to start a plot before running the loop
% Otherwise you'll get "Unrecognized field name "currentAxes"." errors
plot(cS,'faceAlpha',0.5)

t = tiledlayout(8,12,'TileSpacing','tight','Padding','tight',...
    'TileIndexing', 'rowmajor');
for k = 1:length(sS)
  ax = nexttile;
  plot(cS,'faceAlpha',0.5,'parent',ax)
  title(ax,['\textbf{' int2str(k) '}:' char(sS(k).n,'latex')],'Interpreter','latex')
  axis on
  xlabel X; ylabel Y; zlabel Z;
  hold on
  plot(cS,sS(k),'facecolor','red','parent',ax)
  %plottingConvention.default3D().setView
  % Load direction
  % arrow3d(0.4*xvector,'faceColor','red','linewidth',3)
  hold off
end

%% 81 - Plot all the m vectors (with rotation) - ADAM
fprintf('Processing figure (81) -  Plot all the m vectors (with rotation) - ADAM\n');


figure(81);
% figure;

% For some reason, need to start a plot before running the loop
% Otherwise you'll get "Unrecognized field name "currentAxes"." errors
plot(ori*cS,'faceAlpha',0.5);
xlim([-0.5 0.5]); ylim([-0.45 0.45])
t = tiledlayout(8,12,'TileSpacing','tight','Padding','tight',...
    'TileIndexing', 'rowmajor');
for k = 1:length(sS)
  ax = nexttile;
  plot(ori*cS,'faceAlpha',0.5,'parent',ax)
  title(ax,['\textbf{' int2str(k) '}:' char(sS(k).n,'latex')],'Interpreter','latex')
  axis on
  xlabel X; ylabel Y; zlabel Z;
  hold on
  plot(ori*cS,ori*sS(k),'facecolor','red', 'arrowLineWidth', 3,'LineWidth',1.5,'parent',ax)
  %plottingConvention.default3D().setView
  % Load direction
  arrow3d(0.4*xvector,'faceColor','k','linewidth',3)
  hold off
end

%% 811 - Plot all the m vectors (with rotation) - CP
fprintf('Processing figure (81) -  Plot all the m vectors (with rotation) - CP\n');

% Generate a colormap (you can use any MATLAB colormap here)
numVectors = length(sS);
% cmap = jet(numVectors); % Use 'parula' colormap; replace with 'hot', 'cool', etc.
cmap = flipud(jet(numVectors));

figure(811);
% figure;
% For some reason, need to start a plot before running the loop
% Otherwise you'll get "Unrecognized field name "currentAxes"." errors
plot(ori * cS, 'faceAlpha', 0.5);
set(gcf, 'units', 'pixels', 'position', [13,4,1610,837]);
xlim([-0.5 0.5]); ylim([-0.45 0.45]);

t = tiledlayout(8, 12, 'TileSpacing', 'tight', 'Padding', 'tight', ...
    'TileIndexing', 'rowmajor');

% Loop through all vectors in sS
for k = 1:numVectors
    
    ax = nexttile;
    
    % Plot the crystal shape
    plot(ori * cS, 'faceAlpha','facecolor',[0.1216,0.9804 , 0.3098], 0.2, 'parent', ax ,'LineWidth', 1);
    
    % Title for each plot
    title(ax, ['\textbf{' int2str(k) '}:' char(sS(k).n, 'latex')], ...
          'Interpreter', 'latex');
    axis on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    

    hold on;
    
    % % Plot the rotated shape or vector with a color from the colormap
    % plot(ori * cS, ori * sS(k), 'facecolor', cmap(k, :),0.8, ...
    %      'arrowLineWidth', 3, 'LineWidth', 1.5, 'parent', ax);

        % Plot the rotated shape or vector with a color from the colormap
    plot(ori * cS, ori * sS(k), 'facecolor','b', ...
         'arrowLineWidth', 3, 'LineWidth', 1.5, 'parent', ax);
    
    % % Add an arrow (customizable appearance)
    % arrow3d(0.4 * xvector, 'faceColor', 'k', 'linewidth', 3);
    
    hold off;
end

%% 82 - Plot all the m vectors (with rotation and schmid factor) - ADAM

% Assume uniaxial tension in x direction
sigma = stressTensor.uniaxial(xvector);

% rotate the slipSystem to EBSD axis
% Gives a warning about rotating in specimen coordinates
% sS_rot=ori*sS gives same warning...
% But the pole figure rotation looks correct for this order
sS_rot= slipSystem(ori*b_3dvec, ori*m_3dvec);

% take absolute magnitude, otherwise -x stress != x stress
tau_rot=abs(sS_rot.SchmidFactor(sigma))
%tau=ori*sS.SchmidFactor(sigma)

[tauMax,id] = sort(tau_rot,'descend')

figure(82);
%fig = gcf;
%ax = fig.CurrentAxes;
% For some reason, need to start a plot before running the loop
% Otherwise you'll get "Unrecognized field name "currentAxes"." errors
plot(ori*cS,'faceAlpha',0.5)

t = tiledlayout(8,12,'TileSpacing','tight','Padding','tight',...
    'TileIndexing', 'rowmajor');
for k = 1:length(id)
  ax = nexttile;
  plot(ori*cS,'faceAlpha',0.5,'parent',ax)
  title(ax,['\textbf{' int2str(id(k)) '}:' num2str(tauMax(k))],'Interpreter','latex')
  axis on
  xlabel X; ylabel Y; zlabel Z;
  hold on
  plot(ori*cS,ori*sS(id(k)),'facecolor','red','parent',ax)
  %plottingConvention.default3D().setView
  % Load direction
  % arrow3d(0.4*xvector,'faceColor','red','linewidth',3)
  hold off
end
%% 822 - Plot all the m vectors (with rotation and schmid factor) - CP
fprintf('Processing figure (82) -  Plot all the m vectors (with rotationand schmid factor) - CP\n');

% Assume uniaxial tension in x direction
% sigma = stressTensor.dotial(xvector)
sigma = stressTensor.uniaxial(xvector);
cmap = flipud(jet(numVectors));

% rotate the slipSystem to EBSD axis
% Gives a warning about rotating in specimen coordinates
% sS_rot=ori*sS gives same warning...
% But the pole figure rotation looks correct for this order
sS_rot= slipSystem(ori*b_3dvec, ori*m_3dvec);

% take absolute magnitude, otherwise -x stress != x stress
tau_rot=abs(sS_rot.SchmidFactor(sigma))
%tau=ori*sS.SchmidFactor(sigma)

[tauMax,id] = sort(tau_rot,'descend')

figure(822);
%fig = gcf;
%ax = fig.CurrentAxes;
% For some reason, need to start a plot before running the loop
% Otherwise you'll get "Unrecognized field name "currentAxes"." errors
plot(ori*cS,'faceAlpha','facecolor',[0.1216,0.9804 , 0.3098],0.2,'LineWidth', 1.2)
set(gcf, 'units', 'pixels', 'position', [13,4,1610,837]);
% xlim([-0.45 0.45]); ylim([-0.35 0.35])
xlim([-0.5 0.5]); ylim([-0.45 0.45]);
t = tiledlayout(8,12,'TileSpacing','tight','Padding','tight',...
    'TileIndexing', 'rowmajor');
for k = 1:length(id)
  ax = nexttile;
  plot(ori*cS,'faceAlpha','facecolor',[0.1216,0.9804 , 0.3098],0.2,'parent',ax,'LineWidth', 1.2)
  title(ax,['\textbf{' int2str(id(k)) '}:' num2str(tauMax(k))],'Interpreter','latex')
  axis on; xlim([-0.45 0.45]); ylim([-0.35 0.35])

  xlabel X; ylabel Y; zlabel Z;
  hold on
  plot(ori*cS,ori*sS(id(k)),'facecolor', cmap(k, :),1.5,'parent',ax,'LineWidth', 1.5)
  % plot(ori*cS,ori*sS(id(k)),'facecolor', 'k',0.8,'parent',ax,'LineWidth', 1.5)

  %plottingConvention.default3D().setView
  % Load direction
  arrow3d(0.4*xvector,'faceColor','red','linewidth',3)
  hold off
end
axis tight
%% 83 - On pole figure, with markersize a function of work - Adam

figure(82)

scatter(ori*m_3dvec,'grid','on','antipodal',...
    'MarkerSize',30,'Marker','x','MarkerEdgeColor', 'k')
hold on
scatter(ori*sS(id).n,...
    'MarkerSize',tauMax*1000,'grid','on','antipodal')
hold off
%%83 - Pole figure with plane - CP
%%=== Center Plane_with_BothArrows.fig at (0,0) on the PF (figure 83) ===
figure(82); 
axPF = gca; 
hold(axPF,'on');

planePath = '/Users/celesteperez/Desktop/Plane_with_BothArrows.fig';

% -- Open plane fig invisibly
srcFig = openfig(planePath,'new','invisible');
srcAx  = findobj(srcFig,'Type','axes','-depth',1);
set(srcFig,'Color','w');   % use solid bg; we’ll handle alpha ourselves
if ~isempty(srcAx), set(srcAx,'Color','w'); end

% -- Export to PNG (no reliance on built-in transparency)
tmpPng = fullfile(tempdir,'plane_overlay_pf.png');
ok = true;
try
    if ~isempty(srcAx)
        exportgraphics(srcAx,tmpPng,'ContentType','image','Resolution',300);
    else
        exportgraphics(srcFig,tmpPng,'ContentType','image','Resolution',300);
    end
catch
    ok = false;
end
if ~ok
    % Fallback via getframe
    try
        fr = getframe(srcAx);
    catch
        fr = getframe(srcFig);
    end
    [rgb,~] = frame2im(fr);
    imwrite(rgb,tmpPng);
end
close(srcFig);

% -- Read PNG (+ alpha if present)
[img, map, alpha] = imread(tmpPng);
if ~isempty(map), img = ind2rgb(img,map); end   % -> double [0..1]
if ~isa(img,'double'), img = im2double(img); end

% Build AlphaData: prefer embedded alpha; else use uniform 0.82
if ~isempty(alpha)
    A = double(alpha)/255;                      % size: [H W]
else
    A = 0.82;                                   % scalar fallback
end

% If AlphaData must be a matrix, expand scalar to [H W]
[imh, imw, ~] = size(img);
if isscalar(A)
    A = repmat(A, imh, imw);
end

% --- Compute a centered box at (0,0) in PF data units
xl = xlim(axPF); yl = ylim(axPF);
R  = 0.5 * min(diff(xl), diff(yl));   % PF radius
scale = 0.85;                         % 0<scale<=1 (tweak size)
ar = imw/imh;                         % image aspect ratio (w/h)

% Choose width/height to fit inside PF circle, centered at (0,0)
w = 2*R*scale; 
h = w/ar;
if h/2 > R*scale
    h = 2*R*scale; 
    w = h*ar;
end
x1 = -w/2; x2 =  w/2; 
y1 = -h/2; y2 =  h/2;

% --- Draw on the existing PF axes (no new axes)
hImg = image(axPF, 'XData',[x1 x2], 'YData',[y1 y2], 'CData',img);
% Robust AlphaData assignment (avoid size errors)
try
    set(hImg,'AlphaData',A);
catch
    % Last-ditch: use uniform alpha
    set(hImg,'AlphaData',.84);
end
set(hImg,'HitTest','off');   % don’t steal clicks
uistack(hImg,'top');

hold(axPF,'off');


 %% 1000 Adam - Original (orientation is not correct)
    figure(1000);
plot(cS,'faceAlpha',0.5)
hold on
plot(cS,sS(47),'facecolor','blue','label','b')
arrow3d(-0.8*sS(47).n,'faceColor','black','linewidth',1,'label','n')
plottingConvention.default3D().setView

%arrow3d(0.4*r,'faceColor','red','linewidth',2,'label','r')
hold off
 
%% Plot on crystal shape
% https://mtex-toolbox.github.io/CrystalShapes.html 

%% Figure 1001 - b direction is wrong
idx=7; % Grain 4 with prior labelling
% idx=9; % Grain 5 with prior labelling
% idx=3; % Grain 7 with prior labelling
HPVnum = 96;

ori=grains.meanOrientation(idx);
%%Draw individual grains
figure(1001)

plot(cSGrains(idx),'facecolor',[0.1216,0.9804 , 0.3098],0.8,'coordinates','on','faceAlpha',0.5,'LineWidth', 3);
axis on; hold on
xlabel X; ylabel Y; zlabel Z;
plot(cSGrains(idx),sS(HPVnum),'facecolor','k','label','b'); 
% 
arrow3d(-0.8*sS(HPVnum).n,'faceColor','black','linewidth',2,'faceAlpha',0.5,'label','n')

plottingConvention.default3D().setView
hold off

% Set the 3D view (azimuth, elevation)
% (0, 270) gives the X-Y plane with Y positive downward
% need to match what the EBSD plotting axes display
view(0,90);
% view(3)
return
%% 73 - Plot on crystal shape
% https://mtex-toolbox.github.io/CrystalShapes.html 
% https://mtex-toolbox.github.io/SlipSystems.html
ebsd_filled =ebsd ;
figure(73)

 fprintf('Fig (73) - Filled EBSD \n'); 
[~, ebsd_filled.grainId] = calcGrains(ebsd_filled('indexed'), 'angle', 3*degree);
ebsd_filled = fill(ebsd_filled('indexed'), grains);clf;
plot(ebsd_filled(phase_name), ebsd_filled(phase_name).orientations,...
    'figSize', 'large','coordinates','on');
hold on;
% DRAW GRAIN BOUNDARY
plot(grains.boundary, 'linewidth', 2);
hold on

% Choose a grain
grainID = 7;
idx = find(grains.id == grainID);
ori = grains(idx).meanOrientation;
center = grains(idx).centroid;

% Scale crystal
scale = 0.7 * sqrt(grains(idx).area);
cSGrain = ori * cS * scale;

% Transform the slip system to sample coordinates
sS_rot = ori * sS;

% Plot crystal shape on top of EBSD map
plot(center + cSGrain, 'faceColor', [1, 0.6, 0.6], ...
    'edgeColor', 'k', 'LineWidth', 2, 'faceAlpha', 0.7); hold on;

varinum = 7;
% plot(center + cSGrain, sS_rot(varinum), ...
%     'FaceColor', 'blue'); 


arrow3d(vector3d(center, center + 10 * sS_rot.b), ...
    'FaceColor', 'black', 'label', 'b'); hold off;

% 
% % Add arrow for slip plane normal (in sample coordinates)
% arrow3d(vector3d(center, center + 20 * sS_rot(varinum).n), ...
%     'FaceColor', 'magenta', 'label', 'n'); 
%     hold off;

%% Try using schmid factor for uniaxial

sigma = stressTensor.uniaxial(xvector)

tau=sS.SchmidFactor(sigma)

[tauMax,id] = max((tau))

%% in plane angles

% transpose(rad2deg(acos(dot(cross(ori*sS(id).n,zvector),yvector))))
%% Save version information

Version_output("Version_Flag.txt")

%% --------------------- Trying to plot layers on cube 09/2025 ---------------------
%% 81 Plot all the m vectors (with rotation) - CP
warning ('off')

% Generate a colormap (you can use any MATLAB colormap here)
numVectors = length(sS);
cmap = jet(numVectors); % Use 'parula' colormap; replace with 'hot', 'cool', etc.

% figure(81);
figure;
% For some reason, need to start a plot before running the loop
% Otherwise you'll get "Unrecognized field name "currentAxes"." errors
plot(ori * cS, 'faceAlpha', 0.5);
set(gcf, 'units', 'pixels', 'position', [13,4,1610,837]);
xlim([-0.5 0.5]); ylim([-0.45 0.45]);


% Loop through all vectors in sS
for k = 1:numVectors
    ax = nexttile;
    
    % Plot the crystal shape
    plot(ori * cS, 'faceAlpha','facecolor',[0.1216,0.9804 , 0.3098], 0.2, 'parent', ax ,'LineWidth', 1);
    
    % Title for each plot
    title(ax, ['\textbf{' int2str(k) '}:' char(sS(k).n, 'latex')], ...
          'Interpreter', 'latex');
    axis on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    

    hold on;
    
    % % Plot the rotated shape or vector with a color from the colormap
    % plot(ori * cS, ori * sS(k), 'facecolor', cmap(k, :),0.8, ...
    %      'arrowLineWidth', 3, 'LineWidth', 1.5, 'parent', ax);
     plot(ori * cS, ori * sS(k), 'facecolor', 'b',0.8, ...
         'arrowLineWidth', 3, 'LineWidth', 1.5, 'parent', ax);
    % % Add an arrow (customizable appearance)
    % arrow3d(0.4 * xvector, 'faceColor', 'k', 'linewidth', 3);
    
    hold off;
end

%% Plots Rotated Lattice w/ HPV - works
idx=7; % Grain 4 with prior labelling
% idx=9; % Grain 5 with prior labelling
% idx=3; % Grain 7 with prior labelling
% close all
ori=grains.meanOrientation(idx);
%%Draw individual grains

% range = numVectors;
range = [ 51:93 ];
% Loop through all vectors in sS
for n = 1:length(range)
    k = range(n)
    close all
    figure(100+k)
    % figure
    % figure('Name','Variant #',k)

    % % Plot the crystal shape
    plot(ori * cS, 'faceAlpha','facecolor',[0.1216,0.9804 , 0.3098], 0.2, 'LineWidth', 1);

    % Title for each plot
    title(['\textbf{' int2str(k) '}:' char(sS(k).n, 'latex')], ...
          'Interpreter', 'latex', 'fontsize', 20);
    axis on;
    xlabel('X'); ylabel('Y'); zlabel('Z');


    hold on;

    % Plot the rotated shape or vector with a color from the colormap
    plot(ori * cS, ori * sS(k), 'facecolor', 'b',0.8, ...
         'arrowLineWidth', 3, 'LineWidth', 1.5);
     

    % % Add an arrow (customizable appearance)
    arrow3d(0.4 * xvector, 'faceColor', 'r', 'linewidth', 1.5);
    
    hold off;
    view(0,90);
% end
%% 201 - Draw plane only - WORKS
figure(101); clf; hold on
title(['\textbf{' int2str(k) '}:' char(sS(k).n, 'latex')], ...
          'Interpreter', 'latex', 'fontsize', 20);

% 1) Init MTEX axes (keep handle so we don't nuke our own patches later)
hCube = plot(ori * cS, 'faceAlpha', 0, 'edgeAlpha', 0);  % invisible cube

% 2) Arrow (and MTEX’s own minimal plane glyph)
plot(ori * cS, ori * sS(k), ...
     'faceColor', cmap(k,:), 0.8, ...
     'arrowLineWidth', 3, ...
     'arrowFaceColor', 'k', ...
     'arrowEdgeColor', 'k', ...
     'LineWidth', 1.5); hold on
% 2) Arrow (and MTEX’s own minimal plane glyph)
plot(ori * cS, ori * sS(k), ...
     'faceColor','b', 0.8, ...
     'arrowLineWidth', 3, ...
     'arrowFaceColor', 'k', ...
     'arrowEdgeColor', 'k', ...
     'LineWidth', 1.5); hold on;
% % Add an arrow (customizable appearance)
    arrow3d(0.4 * xvector, 'faceColor', 'r', 'linewidth', 1);
    hold off;

axis off
% xlabel('X'); ylabel('Y'); zlabel('Z');
%% Creates .fig of HPV Plane
figure(101); clf; hold on

fig = gcf; ax = gca;

title(ax, ['\textbf{' int2str(k) '}:' char(sS(k).n, 'latex')], ...
      'Interpreter','latex','FontSize',20);

% Ensure the orientation crystal symmetry matches the crystal shape
% (prevents the specimen/crystal coords warning)
if exist('cS','var') && isa(cS,'crystalSymmetry') && isa(ori,'orientation')
    if ~isequal(ori.CS, cS); ori = orientation(ori, cS, ori.SS); end
end


% 1) Invisible cube to initialize MTEX axes
plot(ori*cS, 'faceAlpha',0, 'edgeAlpha',0, 'Parent', ax);

% 2) Draw the plane + MTEX arrow (use proper name-value pairs; add tags)
plot(ori*cS, ori*sS(k), ...
     'faceColor', cmap(k,:), 'faceAlpha', 0.8, ...
     'arrowLineWidth', 3, ...
     'arrowFaceColor', 'k', ...
     'arrowEdgeColor', 'k', ...
     'LineWidth', 1.5, ...
     'Tag','MTEX_Plane');   % tag will be pushed to created graphics objects

% (Optional second tint layer; also correctly name 'faceAlpha')
plot(ori*cS, ori*sS(k), ...
     'faceColor','b', 'faceAlpha', 0.25, ...
     'arrowLineWidth', 3, ...
     'arrowFaceColor', 'k', ...
     'arrowEdgeColor', 'k', ...
     'LineWidth', 1.0, ...
     'Tag','MTEX_Plane_Tint');

% 3) Your custom 3D arrow (keep its handle)
hArrowCustom = arrow3d(0.4 * xvector, 'faceColor','r', 'linewidth',1);
set(hArrowCustom, 'Tag','CustomArrow');

% 4) Clean-up: make background transparent and hide axes
set(ax,'Color','none'); set(fig,'Color','none');
axis(ax,'equal'); axis(ax,'tight'); axis off
ax.XColor='none'; ax.YColor='none'; ax.ZColor='none';

% 5) Remove everything except plane patches and both arrows
%    - Keep: any object tagged 'MTEX_Plane*' or 'CustomArrow'
kids = allchild(ax);
for i = 1:numel(kids)
    obj = kids(i);
    tg  = get(obj,'Tag');
    if ischar(tg) && (startsWith(tg,'MTEX_Plane') || strcmp(tg,'CustomArrow'))
        continue; % keep
    end
    tp = get(obj,'Type');
    switch tp
        case {'text','line','quiver','quivergroup'}
            delete(obj);
        case {'surface','patch'}
            % If it's not one of the plane/arrow objects, NaN it out (safer than delete)
            if ~(ischar(tg) && (startsWith(tg,'MTEX_Plane') || strcmp(tg,'CustomArrow')))
                if isprop(obj,'ZData') && ~isempty(get(obj,'ZData'))
                    Z = get(obj,'ZData'); set(obj,'ZData',nan(size(Z)));
                end
                if isprop(obj,'CData') && ~isempty(get(obj,'CData'))
                    C = get(obj,'CData'); set(obj,'CData',nan(size(C)));
                end
                if isprop(obj,'Vertices') && ~isempty(get(obj,'Vertices'))
                    V = get(obj,'Vertices'); V(:) = nan; set(obj,'Vertices',V);
                end
            end
        otherwise
            % For any other graphic type, try deleting
            try, delete(obj); end
    end
end

% -------------- 6) Save .fig with transparent background
% % outPath = '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data/Plane_with_BothArrows.fig';
% outPath = sprintf(['/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data/HPV_Planes/' ...
%                    'Plane_with_BothArrows_%d.fig'], k);
% 
% try
%     savefig(fig, outPath);
%     fprintf('✅ Saved to:\n%s\n', outPath);
% catch ME
%     warning('⚠️ Could not save figure:\n%s', ME.message);
% end


end
return
%%
%% =====================================================
%  Plot Grain 4 3D Stack + Import + Resize Plane Figure
% ======================================================
hpv_num = [70];
% --- Load your 3D stack (Grain4_3D_stack.mat) ---
pname = '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data';
fname = [pname '/Grain4_3D_stack.mat'];
load(fname, 'S');

figure(444); clf;
set(gcf, 'Color', 'w');
ax = axes('Parent', gcf); hold(ax,'on');
colormap(gcf, S.cmap);

% --- Plot the 3D stack surfaces and patches ---
for k = 1:numel(S.surfaces)
    surf(S.surfaces{k}.X, S.surfaces{k}.Y, S.surfaces{k}.Z, S.surfaces{k}.C, ...
         'EdgeColor', S.surfaces{k}.EdgeColor);
end
for k = 1:numel(S.patches)
    patch('XData',S.patches{k}.X,'YData',S.patches{k}.Y,'ZData',S.patches{k}.Z, ...
          'FaceColor',S.patches{k}.FaceColor,'FaceAlpha',S.patches{k}.FaceAlpha, ...
          'EdgeColor',S.patches{k}.EdgeColor,'LineWidth',S.patches{k}.LineWidth);
end

axis(ax,'vis3d','equal');
set(ax,'CLim',S.clim,'DataAspectRatio',S.dataaspect);
view(ax,S.view);
campos(ax,S.campos); camtarget(ax,S.camtarget); camva(ax,S.camva);
view(ax,0,90); hold(ax,'on');

% %---------------------------------------------------------------------
% % --- LOAD VERIFICATION PLANE - PlanesOnly.fig invisibly ---
% planeFig = openfig('/Users/celesteperez/Desktop/PlanesOnly.fig','invisible');
% planeAx  = findobj(planeFig,'Type','axes');
% % - Copy all children (planes, etc.) ---
% planeChildren = get(planeAx, 'Children');
% % - Shift every object's XData by -100 ---
% xShift = -300;   % negative = move left, positive = move right
% 
% for i = 1:numel(planeChildren)
%     obj = planeChildren(i);
%     % Handle surfaces, patches, or lines
%     if isprop(obj, 'XData')
%         X = get(obj, 'XData');
%         if ~isempty(X)
%             set(obj, 'XData', X + xShift);
%         end
%     end
% end
% % - Copy shifted objects into your main figure axes ---
% copyobj(planeChildren, ax);
% close(planeFig);   % close the temp figure

%---------------------------------------------------------------------
%%=== LOAD HPV PLANE ===
folder   = '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data/HPV_Planes';
fname    = sprintf('Plane_with_BothArrows_%d.fig', hpv_num);   % or '%02d' if zero-padded
planePath = fullfile(folder, fname);

assert(isfile(planePath), 'File not found: %s', planePath);

tmpFig = openfig(planePath, 'new', 'invisible');
tmpAx  = findobj(tmpFig, 'Type', 'axes', '-depth', 1);




% planePath = '/Users/celesteperez/Desktop/Plane_with_BothArrows.fig';
% tmpFig = openfig(planePath, 'new', 'invisible');
% tmpAx  = findobj(tmpFig, 'Type', 'axes', '-depth', 1);

% Collect plane + arrows (tag-based preferred)
objs = findobj(tmpAx, '-regexp', 'Tag', '^(MTEX_Plane|MTEX_Plane_Tint|CustomArrow)$');
if isempty(objs)
    objs = findobj(tmpAx, '-or', 'Type','patch', 'Type','surface', 'Type','line');
end

%  Copy them into your current 3D figure
newObjs = copyobj(objs, ax);
close(tmpFig);

%  Resize the imported plane to match the 6-point plane ===
%  Target plane bounding box from your coordinates
xRange_target = [3640 3710];
yRange_target = [625  720];
zRange_target = [-130   0];
span_target   = [diff(xRange_target), diff(yRange_target), diff(zRange_target)];

%  Collect XYZ from imported objects
allX = []; allY = []; allZ = [];
for i = 1:numel(newObjs)
    o = newObjs(i);
    if isprop(o,'XData') && isprop(o,'YData') && isprop(o,'ZData')
        X = get(o,'XData'); Y = get(o,'YData'); Z = get(o,'ZData');
        if ~isempty(X) && ~isempty(Y) && ~isempty(Z)
            allX = [allX; double(X(:))];
            allY = [allY; double(Y(:))];
            allZ = [allZ; double(Z(:))];
        end
    end
end

%  Bail if no geometry found
if isempty(allX)
    warning('No valid geometry found in imported plane.');
else
    bbox_plane = [min(allX) max(allX);
                  min(allY) max(allY);
                  min(allZ) max(allZ)];
    span_plane = [diff(bbox_plane(1,:)), diff(bbox_plane(2,:)), diff(bbox_plane(3,:))];
    span_plane(span_plane==0) = 1;

    % Uniform scale factor
    ratios = span_target ./ span_plane;
    scaleFactor = median(ratios(isfinite(ratios)));
    if isempty(scaleFactor) || scaleFactor <= 0, scaleFactor = 1; end

    % Centers and translation
    center_plane  = mean(bbox_plane,2)';
    center_target = [mean(xRange_target), mean(yRange_target), mean(zRange_target)];
    translation   = center_target - center_plane;

    % Apply transform
    for i = 1:numel(newObjs)
        o = newObjs(i);
        if isprop(o,'XData') && isprop(o,'YData') && isprop(o,'ZData')
            X = get(o,'XData'); Y = get(o,'YData'); Z = get(o,'ZData');
            Xn = (double(X) - center_plane(1))*scaleFactor + center_target(1);
            Yn = (double(Y) - center_plane(2))*scaleFactor + center_target(2);
            Zn = (double(Z) - center_plane(3))*scaleFactor + center_target(3);
            set(o,'XData',Xn,'YData',Yn,'ZData',Zn);
        end
    end
    fprintf('Imported plane resized by %.3fx and centered at [%.1f, %.1f, %.1f]\n', ...
        scaleFactor, center_target);
end

%  Finalize visualization ===
axis(ax,'vis3d','equal');

set(ax,'YDir','normal');
title(ax, 'Grain 4 Stack + Imported Plane');
xlim([3600,3750]); ylim([610 750])


% view(ax,-86,80);
view(ax,-30,24);


% ------------ Most aligned HPV
m = [0.742, -0.6445, -0.1844];
m_norm = m / norm(m);

xlsxPath = '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/HPV_MeasuredVSCalculated/HPV_MeasuredVSCalculated.xlsx';
sheet = 'Sheet1';

% Read vectors from B7:D102
V = readmatrix(xlsxPath, 'Sheet', sheet, 'Range', 'B7:D102');

% Normalize each row
Vnorm = V ./ vecnorm(V, 2, 2);     % divide each row by its magnitude

% Compute cos(theta) = dot(m,v)/(||m|| ||v||)
cosTheta = Vnorm * m_norm.';       % N×1

% Find which is most aligned
[maxVal, idx] = max(cosTheta);

fprintf('Most aligned vector is row %d (Excel row %d) with cos(theta)=%.4f\n', ...
    idx, idx + 6, maxVal);  % +6 because data starts at Excel row 7
