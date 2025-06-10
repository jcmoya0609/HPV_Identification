
%% Read in stress data


Data=readmatrix("Example_Data/FEA220620-AAC-009-stress.txt")

%% Plot as a reshape
%Data is 400x80

% reshape the data to match
% Also rescale from MPa
sigma_xx=reshape(Data(:,4),[400,80])/1000000
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