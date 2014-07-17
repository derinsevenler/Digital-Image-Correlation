function [validx,validy]=resume_automate_image(grid_x,grid_y,filenamelist,resultsimcorrx,resultsimcorry);

% Code to resume image correlation after error or user intervention
% Programmed by Chris
% Last revision: 12/26/06

% Load necessary files
if exist('grid_x')==0
    load('grid_x.dat')              % file with x position, created by grid_generator.m
end
if exist('grid_y')==0
    load('grid_y.dat')              % file with y position, created by grid_generator.m
end
if exist('filenamelist')==0
    load('filenamelist')            % file with the list of filenames to be processed
end

if exist('resultsimcorrx')==0
    resultsimcorrx = dlmread('resultsimcorrx.txt','\t');              % file with x position, created by grid_generator.m
end
if exist('resultsimcorry')==0
    resultsimcorry = dlmread('resultsimcorry.txt','\t');              % file with y position, created by grid_generator.m
end

[r c]=size(grid_x);
[Imagenum Rasternum]=size(resultsimcorrx);

if r*c<Rasternum
    resultsimcorrx(:,(r*c+1):Rasternum)=[];
    resultsimcorry(:,(r*c+1):Rasternum)=[];
end

resultsimcorrx(Imagenum,:)=[];
resultsimcorry(Imagenum,:)=[];


validx=resultsimcorrx';
validy=resultsimcorry';

save resultsimcorrx.txt resultsimcorrx -ascii -tabs
save resultsimcorry.txt resultsimcorry -ascii -tabs

[validx,validy]=automate_image(grid_x,grid_y,filenamelist,validx,validy);