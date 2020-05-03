%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_contours()
%
% make contours around the kernel density estimated from a
% list of 2d coordinates and save out the resulting plots
%  
% REQUIRED ARGUMENTS:
% 
% files containing coordinates of neurons must be present in 
% ../make_contours/text/ they must be in in csv format with the last 
% N_LANDMARK rows containing registration points that will not be plotted.

% OPTIONAL ARGUMENTS:
%
% xGrid
% yGrid - limits for grid used for kernel density estimate (microns)
%         defaults: xGrid = [0 1000], yGrid = [-400 600];
% xLims
% yLims - limits for plotting hemicord cartoon (microns)
%         defaults: xLims = [0 650],  yLims = [-400 500];
% 
% parameters for kernel density estimation are hardcoded as constants
% at the beginning of this text file
%
% N_BINS = 256, how big is the N_BINS x N_BINS kernel density map
% N_LEVELS = 9, number of contour levels to plot on kernel density estimate
% N_LANDMARKS = 4, number of registration coordinates present
%                  at the end of each data file
%
% RETURN VALUES:
%
% saves scatter and kernel density plots to ../make_contours/plots/
% allData - cell array with the coordinates of each neuron in each dataset
% names -  filename of each text files read in
%
% tamachado, updated 11/2014
% tam2138@columbia.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updates: 
% 1. The distribution of neurons can now be scaled according to the shape
%    of spinal cord shape specified by the user.
% 2. Files of spinal cord shapes need to be put in ../Contours/ and they 
%    must be in csv format and named as "L3" "C4" or "T4".

% Linjing Fang, updated on 05/23/2018
% lfang@salk.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allData, names] = make_contours()
%% define constants and parameters
% stage = input('Which stage are you plotting? P0, P120 or C120?','s');
% if ~(strcmp(stage,'P0') || strcmp(stage, 'P120') || strcmp(stage, 'C120'))
%     error('Please only type in P0, P120 or C120.');
% end

position = input('Which position are you plotting at? C4, T4 or L3?','s');
if ~(strcmp(position,'C4') || strcmp(position, 'T4') || strcmp(position, 'L3'))
    error('Please only type in C4, T4 or L3.');
end

% constants
CORD_COLOR = [.96 .96 .96]; % color for spinal cord outline (.96 = 245/255)
N_LEVELS = 8;               % how many contour levels should we draw?
N_BINS = 256;               % how big of an NxN matrix to use for kde?
N_LANDMARK = 4;             % these are the last N points in each data file

% add export_fig to the path for saving pretty pdfs if necessary
if ~exist('export_fig','file')
    path = which('make_contours');
    ind = strfind(path,'make_contours.m');
    addpath([path(1:ind-1) 'export_fig']);
end

%% Prepare for plotting
% close all;
% get dataset names
[names, paths] = get_dataset_names;
names=regexprep(names,'_','-');
nDatasets = length(paths);

% accumulate all points as we load datasets; return this value
allData = cell(nDatasets,1);

% scale cord and translate cord outline
addpath(fullfile(pwd,'Contours'));
if strcmp(position,'C4')
   filename = 'C4.csv';
   xLims = [0 930];  yLims = [-570 880]; % limits for plotting cords 
elseif strcmp(position,'T4')
   filename = 'T4.csv';
   xLims = [0 670];  yLims = [-430 660]; % limits for plotting cords  
elseif strcmp(position,'L3')
   filename = 'L3.csv';
   xLims = [0 840];  yLims = [-570 600]; % limits for plotting cords
end


xGrid = [xLims(1)-diff(xLims)/4 xLims(2)+diff(xLims)/4]; % limits for kde xgrid
yGrid = [yLims(1)-diff(yLims)/4 yLims(2)+diff(yLims)/4]; % limits for kde ygrid

delimiterIn = ',';
headerlinesIn = 1;
contour = importdata(filename,delimiterIn,headerlinesIn);
x = contour.data(:,1);
y = contour.data(:,2);
x(abs(x)<10)=0;
x(x<0)=0;
xMin = min(x);
xMax = max(x);
yMin = min(y);
yMax = max(y);
scaled_x = (x-xMin)/(xMax-xMin)*diff(xLims)+xLims(1);
scaled_y = (y-yMin)/(yMax-yMin)*diff(yLims)+yLims(1);
cord = struct('x',scaled_x,'y',scaled_y);

% generate three ticks for each axis
xticks = [xLims(1) 0 xLims(2)];
yticks = [yLims(1) 0 yLims(2)];
if min(xLims) >= 0, xticks = [xLims(1) round(mean(xLims)) xLims(2)]; end
if min(yLims) >= 0, yticks = [yLims(1) round(mean(yLims)) yLims(2)]; end

% get bounds for plotting contours
xBounds = linspace(xGrid(1),xGrid(2),N_BINS);
yBounds = linspace(yGrid(1),yGrid(2),N_BINS);

% set up figures for plotting
contourOverlay = figure; set(gcf,'Color','w');
% contourMatrix = figure; set(gcf,'Color','w');

% colors to use for big matrix of all contours
cols = jet(nDatasets);

%% Make all plots
for ii = 1:nDatasets
    scaled_data = [];
    data=xlsread(paths{ii});
    x_data = data(:,1);
    y_data = data(:,2);
    xMin_data = min(x_data);
    xMax_data = max(x_data);
    yMin_data = min(y_data);
    yMax_data = max(y_data);
    
    % Coordinate Adjustment (scale the data for both x and y axis based 
    % on the same scaling facotr getting from x axis's reference points)
    % (scaled_x(end)-scaled_x(end-2))/(xMax_data-xMin_data) is the scaling
    % factor.
    lscale = (xLims(2)-xLims(1))/(650-0); %scale x to [xLim(2),xLim(1)]
    bscale = (0-yLims(1))/(0-(-400)); %scale y to [yLim(1),0]
%     bscale = (yticks(2)-yLims(1))/(0-(-400)); %scale y to [yLim(1),yticks(2)]
%     bscale = (yLims(2)-yLims(1))/(0-(-400)); %scale y to [yLim(1),yLim(2)]

    scaled_data(:,1) = (x_data-0)*lscale + xLims(1);
    scaled_data(:,2) = (y_data-(-400))*bscale + yLims(1);
    
%     scaled_data(:,1) = (x_data-0)*(scaled_x(end)-scaled_x(end-2))/600 + 0;  
%     scaled_data(:,2) = (y_data-(-400))*(max(scaled_y)-min(scaled_y))/(500-(-400)) + min(scaled_y); 

    % scaled_x(end-2) is the xMax based on the side reference point.
    % scaled_data(:,2) = (y_data-(-400))*(max(scaled_y)-min(scaled_y))/(500-(-400)) + (-400); 
    %scaled_data(:,2) = (y_data-yMin_data)*(scaled_x(end)-scaled_x(end-2))/(xMax_data-xMin_data) + scaled_y(end-1); 
    % scaled_y(end-1) is the yMin based on the bottom reference point.
    
    % save data
    allData{ii} = scaled_data;
    
    % plot data points
    figure; set(gcf,'Color','w');
    cax = subplot(1,2,1);
    make_spinal_cord(cax);
    title(names(ii,:));
    plot(scaled_data(1:end-N_LANDMARK,1),scaled_data(1:end-N_LANDMARK,2),'k.');
    fix_axes(cax)
    
    % plot individual contours based on kernel density
    cax = subplot(1,2,2);
    make_spinal_cord(cax);
  
    if size(data,1)>6
        % get kernel density function
        [~, density] = kde2d(scaled_data(1:end-N_LANDMARK,:),N_BINS,...
            [xGrid(1) yGrid(1)],[xGrid(2) yGrid(2)]);

        % make and contours around kernel density function
        plot_contours(cax,density,xBounds,yBounds);
        fix_axes(cax);

        % plot all data points on one plot
        figure(contourOverlay); 
        cax = subplot(1,2,1);
        if ii == 1, make_spinal_cord(cax); end
        title('all genes');
        hold on
        plot(scaled_data(1:end-N_LANDMARK,1),scaled_data(1:end-N_LANDMARK,2),...
           '.','Color','k','MarkerSize',8);
        fix_axes(cax)

        % plot all contours on one plot
        cax = subplot(1,2,2);
        if ii == 1, make_spinal_cord(cax); end
        plot_contours(cax,density,xBounds,yBounds,cols(ii,:));
        fix_axes(cax)
    end
end

% Plot all datapoints on one graph
figure; 
cax = subplot(1,2,1);
make_spinal_cord(cax);
title('all genes');
AData = [];
for i=1:length(allData)
    ss = size(allData{i},1);
    AData = [AData; allData{i}(1:ss-5,:)];
end

plot(AData(:,1),AData(:,2),'.','Color','k','MarkerSize',8);
fix_axes(cax)
    
% plot all contours on one plot
cax = subplot(1,2,2);
make_spinal_cord(cax);
[~, density] = kde2d(AData,N_BINS,[xGrid(1) yGrid(1)],[xGrid(2) yGrid(2)]);
plot_contours(cax,density,xBounds,yBounds);
fix_axes(cax)

% save out summary plot
% export_fig(contourOverlay,'pdf','plots\summary.pdf');

%% Helper functions used during plotting

% plot contours on kernel density
function plot_contours(cax,density,xBounds,yBounds,color)
    % setup the axes
    if nargin < 5, color = 'k'; end
    axes(cax);
    
    % use contourf instead of contour because contourf draws smoother lines
    [~,h] = contourf(xBounds,yBounds,density,N_LEVELS,'Color',color);
    
    % get rid of any colors on contour shapes
    allH = allchild(h);
    if isempty(allH)
        allH=h;
    end    
    set(allH,'FaceColor','none'); 
    
    % hide box around contours
    % set(allH(end),'EdgeColor','none'); 
end

% print spinal cord outline
function make_spinal_cord(cax)  
    % setup the axes
    axes(cax);
    hold on;
    set(cax,'TickDir','out','Layer','top');
    % draw the spinal cord outline
    patch(cord.x(1:end-N_LANDMARK),cord.y(1:end-N_LANDMARK),CORD_COLOR,'EdgeColor','k');
    % patch(cord.x,cord.y,CORD_COLOR,'EdgeColor','k');
end

% standardize the axes we're using once we've plotted everything in them
function fix_axes(cax)
    axes(cax);
    axis image; xlabel('\mum'); ylabel('\mum');
    ylim(yLims); xlim(xLims);
    set(gca,'XTick',xticks,'YTick',yticks);
end

end