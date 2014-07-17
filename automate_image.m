function [validx,validy]=automate_image

clear;
close all;

% Initialization. Load necessary files and unpack variables.
load('grid_x.dat')              % file with x position, created by grid_generator.m
load('grid_y.dat')              % file with y position, created by grid_generator.m
% load('grid_a.dat')              % file with initial rotation (all 0s), created by grid_generator.m
load('filenamelist')            % file with the list of filenames to be processed

grid_a = zeros(size(grid_x));
base_points_in = [grid_x grid_y];
input_points_in = [base_points_in grid_a];

num_cp = size(grid_x,1);
numImages=size(filenamelist,1);

validx = zeros(num_cp,numImages-1);
validy = validx;
valida = validx;

% Create figure for visualization and waitbar and start timer

figure;
imshow(filenamelist(1,:))           % show the first image
title('Initial grid for Image Corelation')        % put a title
hold on;
plot(grid_x,grid_y,'g+')            % plot the grid onto the image
hold off;
g = waitbar(0,sprintf('Progress'));
set(g,'Position',[275,50,275,50]);

% Perform DIC
tic % start timer
base = imread(filenamelist(1,:));
for i= 1:(numImages-1)
    input = imread(filenamelist((i+1),:));
    last_input = imread(filenamelist((i),:));
    
    xyainput = cpcorra(input_points_in, base_points_in, input, base, last_input,false);
    validx(:,i) = xyainput(:,1);
    validy(:,i) = xyainput(:,2);
    valida(:,i) = xyainput(:,3);
    input_points_in = xyainput;
    
    % Update image
    imshow(filenamelist(i+1,:));
    hold on;
    plot(grid_x,grid_y,'g+');
    plot(input_points_in(:,1), input_points_in(:,2),'r+');
    hold off;
    
    % Update waitbar
    waitbar(i/(numImages-1),g)
end

% Get Rotation Information from cpcorra
% xyainput = cpcorra(xyainput, base_points_in, imread(filenamelist(end,:)), base, base,true);
% validx(:,end) = xyainput(:,1);
% validy(:,end) = xyainput(:,2);
% valida(:,end) = xyainput(:,3);

close(g)

% save
save validx.dat validx -ascii -tabs
save validy.dat validy -ascii -tabs
save valida.dat valida -ascii -tabs

end
