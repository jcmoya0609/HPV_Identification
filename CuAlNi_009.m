%% Built from of CuAlNi_AC_CP_s5_b1.m
% Commit 67f8be50fc22ec2a816eb0a5f8a1f586fef21f79
% Branch "EBSD-import"

% Plotting wireframe crystal orientation
% Plotting Austenite-Twinned Marteniste (m) vectors
% Plotting rotated m vectors
% Using slip & schmid functions in mtex to do avail work


%% CLEAR ALL PRIOR OUTPUT
clc; clear; close all;



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

%% Fig (1) - PLOT BAND CONTRAST 

ebsd.how2plot.east = xvector;
ebsd.how2plot.north = yvector;
ebsd.how2plot
%%

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
fprintf('EBSD Subplot\n'); figure(14)
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

% PLOT EBSD UNFILTERED WITH GB 
figure(23);

cSGrains = grains.meanOrientation * cS * 0.7* sqrt(grains.area);


color_x = oM_x.orientation2color(ebsd(phase_name).orientations);
plot(ebsd(phase_name), color_x, 'figSize', 'huge','coordinates','on');
hold on;
% PLOT ORIENTATION ON TOP
%plot(ebsd(phase_name), ebsd(phase_name).orientations,...
%    'figSize', 'large');
% PLOT GRAIN BOUNDARIES FROM FILTERED DATA
plot(grains.boundary, 'linewidth', 2);
text(grains,grains.id, 'FontSize',20);
% plot Crystal Shape
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

%% Choose a grain
%idx=7; % Grain 4 with prior labelling
%idx=9; % Grain 5 with prior labelling
idx=3; % Grain 7 with prior labelling

ori=grains.meanOrientation(idx);

%% Draw individual grains
figure(51)

plot(cSGrains(idx),'coordinates','on','faceAlpha',0.2);
axis on; hold on
xlabel X; ylabel Y; zlabel Z;
hold off
% Set the 3D view (azimuth, elevation)
% (0, 270) gives the X-Y plane with Y positive downward
% need to match what the EBSD plotting axes display
view(0,90);


%% Plot m vectors rotated

%WITH ROTATION 
figure(61)
scatter(grains.meanOrientation(idx)*m_3dvec,'grid','on','antipodal');


%% Define using slipSystem in mtex

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

%% Plot all the m vectors (with rotation)

figure(81);
% For some reason, need to start a plot before running the loop
% Otherwise you'll get "Unrecognized field name "currentAxes"." errors
plot(ori*cS,'faceAlpha',0.5);

t = tiledlayout(8,12,'TileSpacing','tight','Padding','tight',...
    'TileIndexing', 'rowmajor');
for k = 1:length(sS)
  ax = nexttile;
  plot(ori*cS,'faceAlpha',0.5,'parent',ax)
  title(ax,['\textbf{' int2str(k) '}:' char(sS(k).n,'latex')],'Interpreter','latex')
  axis on
  xlabel X; ylabel Y; zlabel Z;
  hold on
  plot(ori*cS,ori*sS(k),'facecolor','red','parent',ax)
  %plottingConvention.default3D().setView
  % Load direction
  % arrow3d(0.4*xvector,'faceColor','red','linewidth',3)
  hold off
end

%% 
%% Plot all the m vectors (with rotation and schmid factor)

% Assume uniaxial tension in x direction
sigma = stressTensor.uniaxial(xvector)

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

%% On pole figure, with markersize a function of work 
figure(83)

scatter(ori*m_3dvec,'grid','on','antipodal',...
    'MarkerSize',30,'Marker','x','MarkerEdgeColor', 'k')
hold on
scatter(ori*sS(id).n,...
    'MarkerSize',tauMax*1000,'grid','on','antipodal')
hold off


 %%
figure(1000);

plot(cS,'faceAlpha',0.5)
hold on
plot(cS,sS(1),'facecolor','blue','label','b')
arrow3d(-0.8*sS(1).n,'faceColor','black','linewidth',2,'label','n')
plottingConvention.default3D().setView

%arrow3d(0.4*r,'faceColor','red','linewidth',2,'label','r')
hold off

%% Plot on crystal shape
% https://mtex-toolbox.github.io/CrystalShapes.html 
% https://mtex-toolbox.github.io/SlipSystems.html

figure(73)

 fprintf('Fig (31) - Filled EBSD \n'); 
[~, ebsd_filled.grainId] = calcGrains(ebsd_filled('indexed'), 'angle', 3*degree);
ebsd_filled = fill(ebsd_filled('indexed'), grains);clf;
plot(ebsd_filled(phase_name), ebsd_filled(phase_name).orientations,...
    'figSize', 'large','coordinates','on');
hold on;
% DRAW GRAIN BOUNDARY
plot(grains.boundary, 'linewidth', 2);
hold on

% Choose a grain
grainID = 4;
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
plot(center + cSGrain, sS_rot(varinum), ...
    'FaceColor', 'blue'); 


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

%transpose(rad2deg(acos(dot(cross(ori*sS(id).n,zvector),yvector))))
%% Save version information

Version_output("Version_Flag.txt")





