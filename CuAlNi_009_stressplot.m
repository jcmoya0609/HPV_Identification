
%% Read in stress data
clc; clear;close all

Data_STRESS=readmatrix("Example_Data/FEA220620-AAC-009-stress.txt")

%% Plot as a reshape
%Data is 400x80

% reshape the data to match
% Also rescale from MPa
sigma_xx=reshape(Data_STRESS(:,4),[400,80])/1000000
figure(991)

image(imrotate(sigma_xx,90),'CDataMapping','scaled')
colormap("hot")

% Explicitly set the colormap range
clim([50 150])
% Add the colorbar
colorbar

% Reverse the Y axis to match EBSD
set(gca,'YDir','normal')

% Set an equal scale for pixels
daspect([1 1 1])


%% Subplot
addpath('/Users/celesteperez/Desktop/BUCSEK_LAB_MATLAB/Research_Github/Janice_HPV_Github/Example_Data');


Data_STRESS=readmatrix("Example_Data/FEA220620-AAC-009-stress.txt");
Data_STRAIN=readmatrix("Example_Data/FEA220620-AAC-009-strain.txt");
%%STRESS - Plot as a reshape
%Data is 400x80
% reshape the data to match
% Also rescale from MPa (1000000 factor)
% And by applied stress (Target MPa from load step / 100 MPa in model )
sigma_xx=reshape(Data_STRESS(:,4),[400,80])*(46/100)/(1000000);
figure(991)
set(gcf, 'Color', 'w'); % Set figure background to black
image(sigma_xx,'CDataMapping','scaled')
colormap("hot")
% Explicitly set the colormap range
clim([25 75])
% ---> Colorbar
    cbar = colorbar;
    % Adjust colorbar ticks and appearance
    newTicks = linspace(25, 75, 5); % Example: 5 ticks
    cbar.Ticks = newTicks;
    cbar.FontSize = 32; % Adjust the font size of the ticks
    cbar.Color = 'k'; % Set the tick color to white for visibility
    % Customize colorbar label
    cbarLabel = cbar.Label;
   cbar.Label.String = '\sigma \{220\}';
    cbar.TickLabels = arrayfun(@(x) sprintf('%.0f', x),...
        newTicks, 'UniformOutput', false);
    cbarLabel.FontSize = 32; 
    cbarLabel.FontWeight = 'bold'; 
    cbarLabel.Color = 'k';
    % Adjust axes properties
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'k', 'FontSize', 12); % Black background, white axes
    set(gca, 'XTickLabel', []); 
    set(gca, 'YTickLabel', []);
% Reverse the Y axis to match EBSD
set(gca,'YDir','normal')
ylim([150 350])
% Set an equal scale for pixels
daspect([1 1 1])
% Create a new figure for the subplot
figure('Name', 'Stress and Strain Comparison', 'Color', 'w');

% Define a 1x2 tiled layout
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, 'Finite Element Analysis - 46 MPa', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k');

%%Plot 1: Stress
ax1 = nexttile;
sigma_xx = reshape(Data_STRESS(:,4), [400, 80]) * (46 / 100) / 1000000; % Convert Stress to MPa
imagesc(ax1, sigma_xx, 'CDataMapping', 'scaled');
colormap(ax1, "jet");
clim(ax1, [25 75]); % Set color limits
title(ax1, '\sigma_{220}', 'FontSize', 16, 'FontWeight', 'bold');
set(ax1, 'YDir', 'normal', 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 12);
ylim(ax1, [150 350]);
daspect(ax1, [1 1 1]); % Equal aspect ratio

% Add colorbar for Stress
cbar1 = colorbar(ax1);
newTicks = linspace(25, 75, 5);
cbar1.Ticks = newTicks;
cbar1.TickLabels = arrayfun(@(x) sprintf('%.0f', x), newTicks, 'UniformOutput', false);
cbar1.Label.String = '\sigma_{220} (MPa)';
cbar1.Label.FontSize = 16;
cbar1.Label.FontWeight = 'bold';
cbar1.Label.Color = 'k';
cbar1.FontSize = 12;
cbar1.Color = 'k';
hold on
%%Plot 2: Strain
ax2 = nexttile;
epsilon_xx = reshape(Data_STRAIN(:,4), [400, 80]) * (1000000 * (46 / 100)); % Convert Strain to microstrain
imagesc(ax2, epsilon_xx, 'CDataMapping', 'scaled');
colormap(ax2, "jet");
clim(ax2, [145 550]); % Set color limits
title(ax2, '\mu\epsilon_{220}', 'FontSize', 16, 'FontWeight', 'bold');
set(ax2, 'YDir', 'normal', 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 12);
ylim(ax2, [150 350]);
daspect(ax2, [1 1 1]); % Equal aspect ratio

% Add colorbar for Strain
cbar2 = colorbar(ax2);
newTicks = linspace(145, 550, 5);
cbar2.Ticks = newTicks;
cbar2.TickLabels = arrayfun(@(x) sprintf('%.0f', x), newTicks, 'UniformOutput', false);
cbar2.Label.String = '\mu\epsilon_{220}';
cbar2.Label.FontSize = 16;
cbar2.Label.FontWeight = 'bold';
cbar2.Label.Color = 'k';
cbar2.FontSize = 12;
cbar2.Color = 'k';