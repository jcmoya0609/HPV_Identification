%% CLEAR ALL PRIOR OUTPUT
clc; clear; close all;

%% SET-UP PATH 
% only needed initialy

%setMTEXpref('voronoiMethod','jcvoronoi');
%mtexPath= '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/mtex-6.2.beta.3';  % Path to mtex folder
%addpath(mtexPath);  startup_mtex
%addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB');

%% CRYSTAL AND SPECIMEN AND SYMMETRIES 

% CRYSTAL SYMMETRY
CS = {'notIndexed',crystalSymmetry('m-3m', [5.8 5.8 5.8], ...
    'mineral', 'CuAlNi-beta', 'color', [0.53 0.81 0.98])};

% PLOTTING CONVENTION
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','IntoPlane');

%% SPECIFY FILE NAMES 

%RELATIVE PATH TO FILES
pname = './Example_Data';
fname = [pname '/CuAlNi_AC_s5_b1 Specimen 1 Site 1 Map Data 1.h5oina'];
% SAVE NAME AND PATH
plotname='CuAlNi_AC_s5_b1';
savepath='CuAlNi_AC_s5_b1_plots';
phase_name='CuAlNi-beta';

%% IMPORT THE DATA 

ebsd = EBSD.load(fname,CS,'interface','h5oina');
ebsd = rotate(ebsd,rotation.byAxisAngle(xvector,180*degree))

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

fprintf('Processing figure (1) - Band Contrast\n');
figure(1);
%plot(ebsd,ebsd.bc,'micronbar','off');
plot(ebsd,ebsd.bc);
colormap gray; % this makes the image grayscale
mtexColorbar;
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
cS = crystalShape.cube(ebsd.CS)
plot(cS,'faceAlpha',0.2)
drawNow(gcm,'final')
axis on

%% Fig (11, 12, 13, 14) - IPF MAPS

%%================ IPF X ================= %%
fprintf('Processing figure (11) - EBSDX\n');
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

%%================ IPF Maps Subplot ================= %%
fprintf('EBSD Subplot\n'); figure(14)
mtexFig = newMtexFigure('layout', [2, 2]);
phase_name_x = phase_name;  phase_name_y = phase_name;  phase_name_z = phase_name;
% IPF-X
oM_x = ipfHSVKey(ebsd('indexed'));  oM_x.inversePoleFigureDirection = xvector;
color_x = oM_x.orientation2color(ebsd(phase_name_x).orientations);
plot(ebsd(phase_name_x), color_x);  title('IPF-X');  nextAxis;
% IPF-Y
oM_y = ipfHSVKey(ebsd('indexed'));  oM_y.inversePoleFigureDirection = yvector;
color_y = oM_y.orientation2color(ebsd(phase_name_y).orientations);
plot(ebsd(phase_name_y), color_y);  title('IPF-Y');  nextAxis;
% IPF-Z
oM_z = ipfHSVKey(ebsd('indexed'));  oM_z.inversePoleFigureDirection = zvector;
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

ebsd_orig = ebsd; % Grain Reconstruction
ebsd_unfiltered = ebsd; % Start from original
ebsd_filled = ebsd;
ebsd_roi = ebsd; 
ebsd_color = ebsd;

% COLOR FROM ORIGINAL DATA UNFILTERED
color_unfiltered = oM.orientation2color(ebsd_color(phase_name).orientations);

%% GRAIN RECONSTRUCTION

% INITIAL GRAIN RECONSTRUCTION
[grains, ebsd_orig.grainId] = calcGrains(ebsd_orig('indexed'), 'angle', 3*degree);
% REMOVE SMALL GRAINS
ebsd_orig(grains(grains.numPixel < 1000)) = [];
% RECOMPUTE GRAINS AFTER REMOVAL
[grains, ebsd_orig.grainId] = calcGrains(ebsd_orig('indexed'), 'angle', 3*degree);
% SMOOTH THE GRAIN BOUNDARIES (OPTIONAL)
grains = smooth(grains, 5);

%% Di
figure(21);
% Display average orientation of remaining grains
disp(grains.meanOrientation);

%% Fig (22) -  EBSD + GRAIN BOUNDARIES (CP)

% PLOT EBSD UNFILTERED WITH GB 
figure(22); clf;
plot(ebsd_unfiltered(phase_name), color_unfiltered, 'figSize', 'huge');
hold on;
% PLOT ORIENTATION ON TOP
plot(ebsd_unfiltered(phase_name), ebsd_unfiltered(phase_name).orientations,...
    'figSize', 'large');
% PLOT GRAIN BOUNDARIES FROM FILTERED DATA
plot(grains.boundary, 'linewidth', 2);
hold off;

%%=============  Save Figure if Needed =============== %%
if save_files == true
    saveFigure(strcat(savepath, plotname, '-IPF-EBSDZ-3degGB-Raw.png'));
end

%% SELECTING ROI AREA AND CRYSTAL SHAPE

% SELECT ROI
xmin = 570; ymin = -200; dx = 100; dy = 100;
% APPLY ROI
region = [xmin ymin dx dy];
fprintf('Processing figure (31) - Fill EBSD\n');
% CRYSTAL SHAPE
isBig = grains.numPixel>50;
cSGrains = grains(isBig).meanOrientation * cS * 0.7 * sqrt(grains(isBig).area);
%FROM JANICE'S CODE:
%CTM = readmatrix("Shield_CTM_CuAlNi_Results.csv");
%HPVNum = CTM(:,1);
%b = CTM(:,4:6);
%m = CTM(:,7:9); 
%m_3dvec = vector3d(m');
% b_3dvec = vector3d(b');
% m_3dvec=reshape(double(m_3dvec),[number_of_interfaces,3]); %csv
% b_3dvec=reshape(double(b_3dvec),[number_of_interfaces,3]); %csv 
Selected_var = m_3dvec(1);%CTM(1,7:9)

%% Fig (31) - FILLED EBSD + GRAIN BOUNDARIES (CP)
oM_top = ipfHSVKey(CS{2});
oM.inversePoleFigureDirection = zvector;

color_unfiltered = oM.orientation2color(ebsd_color(phase_name).orientations);
% color_unfiltered = oM.orientation2color(grains.meanOrientation);

% PLOT FILLED EBSD
figure(31); fprintf('Fig (31) - Filled EBSD \n'); 
[~, ebsd_filled.grainId] = calcGrains(ebsd_filled('indexed'), 'angle', 3*degree);
ebsd_filled = fill(ebsd_filled('indexed'), grains);clf;
plot(ebsd_filled(phase_name), ebsd_filled(phase_name).orientations, 'figSize', 'large');
hold on;
% DRAW RECTANGLE
rectangle('position',region,'edgecolor','r','linewidth',2); hold on;
% DRAW GRAIN SHAPE
plot(grains(isBig).centroid + cSGrains,'FaceColor',...
    color_unfiltered(isBig,:),'linewidth',2,'FaceAlpha',0.7); hold on;
% DRAW 3D Vector
Z = X.*exp(-X.^2 - Y.^2);
[U,V,W] = surfnorm(X,Y,Z);
quiver3(Selected_var)


% DRAW GRAIN BOUNDARY
plot(grains.boundary, 'linewidth', 2);
hold off;
drawNow(gcm,'final')

% LIST GRAIN ORIENTATION
disp(grains.meanOrientation);

% SAVE IF ENABLED
if save_files == true
    saveFigure(strcat(savepath, plotname, '-IPF-EBSDZ-3degGB-Filled.png'));
end

%% Fig 32 - POLE FIGURE WITH GRAINS
% Not working well, maybe a for loop is needed to plot as PDF
%figure (32);

%color_subset = oM.orientation2color(grains(isBig).meanOrientation);
%oM.orientation2color

%title('Pole Figure - CuAlNi-beta');
%plotPDF(ori, h, 'antipodal', 'MarkerSize', 4); 
%plot(grains(isBig).meanOrientation,grains(isBig).meanOrientation * cS * 0.3,...
%    'FaceColor',color_subset,'add2all') 


 %% Fig (41) -  PLOT ROI 
% 
% % APPLY ROI
% xmax = xmin + dx; ymax = ymin + dy;
% x = [xmin xmax xmax xmin];
% y = [ymin ymin ymax ymax];
% polygon = [x' y'];
% ind = inpolygon(ebsd, polygon);
% ebsd_roi = ebsd_roi(ind);
% 
% % ROI GRAIN RECONSTRUCTION
% [grains, ebsd_roi.grainId] = calcGrains(ebsd_roi('indexed'), 'angle', 3*degree);
% ebsd_roi(grains(grains.grainSize < 1000)) = [];
% [grains, ebsd_roi.grainId] = calcGrains(ebsd_roi('indexed'), 'angle', 3*degree);
% grains = smooth(grains, 5);
% 
% % PLOR ROI EBSD
% figure(41); fprintf('Processing figure (41) - Fill EBSD\n');
% [~, ebsd_roi.grainId] = calcGrains(ebsd_roi('indexed'), 'angle', 3*degree);
% ebsd_filled = fill(ebsd_roi('indexed'), grains);clf;
% plot(ebsd_filled(phase_name), ebsd_filled(phase_name).orientations, 'figSize', 'large');
% hold on;
% plot(grains.boundary, 'linewidth', 2);
% hold off;
% 
% % SAVE IF ENABLED
% if save_files == true
%     saveFigure(strcat(savepath, plotname, '-IPF-EBSDZ-3degGB-Filled.png'));
% end


%% Plot m vectors

figure(51)
scatter(m_3dvec,'grid','on','antipodal')
figure(52)
scatter(m_3dvec,'grid','on')

figure(53)
scatter(m_3dvec(1),'grid','on','antipodal');
annotate(vector3d(1,0,0),'label',{'X'});
annotate(vector3d(0,1,0),'label',{'Y'});
annotate(vector3d(0,0,1),'label',{'Z'});

%% Plot m vectors rotated


figure(61)
scatter( grains.meanOrientation(4)*m_3dvec,'grid','on','antipodal')
figure(62)
scatter( grains.meanOrientation(4)*m_3dvec,'grid','on')

%% Plot on crystal shape
% https://mtex-toolbox.github.io/CrystalShapes.html 
% https://mtex-toolbox.github.io/SlipSystems.html

sS=slipSystem(b_3dvec,m_3dvec)

%sS = [slipSystem.pyramidal2CA(ebsd.CS), ...
%  slipSystem.pyramidalA(ebsd.CS)]
figure(71)
plot(cS,'faceAlpha',0.2)
axis on;
xlabel X;
ylabel Y;
zlabel Z;
hold on
plot(cS,sS(1),'faceColor','blue')
hold off
drawNow(gcm,'final')


%% Save version information

Version_output("Version_Flag.txt")





