%% CLEAR ALL PRIOR OUTPUT
clc; clear; close all;

%% SET-UP PATH 
% only needed initialy

setMTEXpref('voronoiMethod','jcvoronoi');
mtexPath= '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/mtex-6.2.beta.3';  % Path to mtex folder
addpath(mtexPath);  startup_mtex
addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data');

%% CRYSTAL AND SPECIMEN AND SYMMETRIES 

% CRYSTAL SYMMETRY
CS = {'notIndexed',crystalSymmetry('m-3m', [5.8 5.8 5.8], ...
    'mineral', 'CuAlNi-beta', 'color', [0.53 0.81 0.98])};

% PLOTTING CONVENTION
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','IntoPlane');

%% SPECIFY FILE NAMES 

%RELATIVE PATH TO FILES

% pname = './Example_Data'; % Adam
pname = '/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data'; % Celeste
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
plot(cS,'faceAlpha',0.2); hold on; 

%% GRAIN RECONSTRUCTION
ebsd_orig = ebsd;
ebsd_color = ebsd;
ebsd_filled = ebsd;
% INITIAL GRAIN RECONSTRUCTION
[grains, ebsd_orig.grainId] = calcGrains(ebsd_orig('indexed'), 'angle', 3*degree);
% REMOVE SMALL GRAINS
ebsd_orig(grains(grains.numPixel < 1000)) = [];
% RECOMPUTE GRAINS AFTER REMOVAL
[grains, ebsd_orig.grainId] = calcGrains(ebsd_orig('indexed'), 'angle', 3*degree);
% SMOOTH THE GRAIN BOUNDARIES (OPTIONAL)
grains = smooth(grains, 5);

%% CRYSTAL SHAPE
isBig = grains.numPixel>50;
cSGrains = grains.meanOrientation * cS * 0.7 * sqrt(grains.area);

%% Fig (31) - FILLED EBSD + GRAIN BOUNDARIES (CP)
oM_top = ipfHSVKey(CS{2});
oM.inversePoleFigureDirection = zvector;
isBig = grains.numPixel>50;
bigGrainIDs = find(isBig);

color_unfiltered = oM.orientation2color(ebsd_color(phase_name).orientations);
% color_unfiltered = oM.orientation2color(grains.meanOrientation);

% PLOT FILLED EBSD
figure(31); fprintf('Fig (31) - Filled EBSD \n'); 
[~, ebsd_filled.grainId] = calcGrains(ebsd_filled('indexed'), 'angle', 3*degree);
ebsd_filled = fill(ebsd_filled('indexed'), grains);clf;
plot(ebsd_filled(phase_name), ebsd_filled(phase_name).orientations, 'figSize', 'large');
hold on;
% DRAW RECTANGLE
% rectangle('position',region,'edgecolor','r','linewidth',2); hold on;
% DRAW GRAIN SHAPE
plot(grains(isBig).centroid + cSGrains,'FaceColor',...
    color_unfiltered(isBig,:),'linewidth',2,'FaceAlpha',0.7); 
hold on;
% plotInnerFace( cSGrains(4),N,'faceColor','blue')
hold on
% DRAW GRAIN BOUNDARY
plot(grains.boundary, 'linewidth', 2);
hold on
drawNow(gcm,'final')
hold off;

% LIST GRAIN ORIENTATION
disp(grains.meanOrientation);

% SAVE IF ENABLED
if save_files == true
    saveFigure(strcat(savepath, plotname, '-IPF-EBSDZ-3degGB-Filled.png'));
end


%% Plot m vectors rotated

%WITH ROTATION 
figure(61)
scatter( grains.meanOrientation(4)*m_3dvec(7),'grid','on','antipodal')
figure(62)
scatter( grains.meanOrientation(4)*m_3dvec(7),'grid','on')


%% New trial to plot crystal shape and ebsd
% plot an EBSD map
clf % clear current figure
g_x = 500;
g_y = -400;
scaling = 0.7 * sqrt(grains(4).area); % scale the crystal shape to have a nice size
N = Miller(m_3dvec(7),'hkl', CS{2})

figure(41)

% plot at position (500,500) the orientation of the corresponding crystal
plot(ebsd(g_x,g_y).orientations * cS *scaling,'faceAlpha',0.5,'colored','linewidth',2)
hold on
plotInnerFace(ebsd(g_x,g_y).orientations * cS* scaling,N,'faceColor','b')
hold off
axis equal
axis tight

drawNow(gcm,'final')
% h = Miller({1,0,0},{1,1,0},{1,1,1}, CS{2});

%% Fig (31) - FILLED EBSD + GRAIN BOUNDARIES (CP)
oM_top = ipfHSVKey(CS{2});
oM.inversePoleFigureDirection = zvector;
isBig = grains.numPixel>50;
bigGrainIDs = find(isBig);

color_unfiltered = oM.orientation2color(ebsd_color(phase_name).orientations);
% color_unfiltered = oM.orientation2color(grains.meanOrientation);

% PLOT FILLED EBSD
figure(31); fprintf('Fig (31) - Filled EBSD \n'); 
[~, ebsd_filled.grainId] = calcGrains(ebsd_filled('indexed'), 'angle', 3*degree);
ebsd_filled = fill(ebsd_filled('indexed'), grains);clf;
plot(ebsd_filled(phase_name), ebsd_filled(phase_name).orientations, 'figSize', 'large');
hold on;
% DRAW RECTANGLE
% rectangle('position',region,'edgecolor','r','linewidth',2); hold on;
% DRAW GRAIN SHAPE
plot(grains(isBig).centroid + cSGrains,'FaceColor',...
    color_unfiltered(isBig,:),'linewidth',2,'FaceAlpha',0.7); 
hold on;
% plotInnerFace( cSGrains(4),N,'faceColor','blue')
hold on
% DRAW GRAIN BOUNDARY
plot(grains.boundary, 'linewidth', 2);
hold on
drawNow(gcm,'final')
hold off;

% LIST GRAIN ORIENTATION
disp(grains.meanOrientation);

% SAVE IF ENABLED
if save_files == true
    saveFigure(strcat(savepath, plotname, '-IPF-EBSDZ-3degGB-Filled.png'));
end
