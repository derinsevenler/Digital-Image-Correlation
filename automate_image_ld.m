function [validx,validy]=automate_image_ld(displx, disply, grid_x,grid_y,filenamelist,validx,validy);

% Code to start actual image correlation for large displacements, needs the
% displacementin x and y direction
% Programmed by Chris and Rob
% Last revision: 09/10/08

% The automation function is the central function and processes all markers and 
% images by the use of the matlab function cpcorr.m. 
% Therefore the Current directory in matlab has to be the folder where 
%  automate_image.m finds the filenamelist.mat, grid_x.dat and grid_y.dat as well 
% as the images specified in filenamelist.mat. Just type automate_image; and 
% press ENTER at the command line of matlab. 
% At first, automate_image.m will open the first image in the filenamelist.mat and 
% plot the grid as green crosses on top. The next step will need some time since 
% all markers in that image have to be processed for the first image. After correlating 
% image one and two the new raster positions will be plotted as red crosses. On top 
% of the image and the green crosses. The next dialog will ask you if you want to 
% continue with this correlation or cancel. If you press continue, automate_image.m 
% will process all images in the filenamelist.mat. The time it will take to process 
% all images will be plotted on the figure but can easily be estimated by knowing the 
% raster point processing speed (see processing speed). 
% Depending on the number of images and markers you are tracking, this process 
% can take between seconds and days. For 100 images and 200 markers a decent 
% computer should need 200 seconds. To get a better resolution you can always 
% run jobs overnight (e.g. 6000 markers in 1000 images) with higher resolutions. 
% Keep in mind that CORRSIZE which you changed in cpcorr.m will limit your 
% resolution. If you chose to use the 15 pixel as suggested a marker distance of 
% 30 pixel will lead to a full cover of the strain field. Choosing smaller marker 
% distances will lead to an interpolation since two neighboring markers share 
% pixels. Nevertheless a higher marker density can reduce the noise of the strain field.
% When all images are processed, automate_image will write the files validx.mat, 
% validy.mat, validx.txt and validy.txt. The text files are meant to store the result in a 
% format which can be accessed by other programs also in the future.

 
% exist('grid_x')
% exist('grid_y')
% exist('filenamelist')
% exist('validx')
% exist('validy')

% Load necessary files
if exist('displx')==0
    return
end
if exist('disply')==0
    return
end
if exist('grid_x')==0
    load('grid_x.dat')              % file with x position, created by grid_generator.m
end
if exist('grid_y')==0
    load('grid_y.dat')              % file with y position, created by grid_generator.m
end
if exist('filenamelist')==0
    load('filenamelist')            % file with the list of filenames to be processed
end
resume=0;
if exist('validx')==1
    if exist('validy')==1
        resume=1;
        [Rasternum Imagenum]=size(validx);
    end
end


% Initialize variables
input_points_x=grid_x;
base_points_x=grid_x;

input_points_y=grid_y;
base_points_y=grid_y;

if resume==1
    input_points_x=validx(:,Imagenum);
    input_points_y=validy(:,Imagenum);
    inputpoints=1;
end

[row,col]=size(base_points_x);      % this will determine the number of rasterpoints we have to run through
[r,c]=size(filenamelist);                   % this will determine the number of images we have to loop through

% Open new figure so previous ones (if open) are not overwritten
h=figure;
imshow(filenamelist(1,:))           % show the first image
title('Initial Grid For Image Correlation (Note green crosses)')        % put a title
hold on
plot(grid_x,grid_y,'g+')            % plot the grid onto the image
hold off

% Start image correlation using cpcorr.m
g = waitbar(0,sprintf('Processing images'));        % initialize the waitbar
set(g,'Position',[275,50,275,50])                               % set the position of the waitbar [left bottom width height]
firstimage=1;

if resume==1
    firstimage=Imagenum+1
end
    base = uint8(mean(double(imread(filenamelist(1,:))),3));            % read in the base image ( which is always  image number one. You might want to change that to improve correlation results in case the light conditions are changing during the experiment
for i=firstimage:(r-1)               % run through all images
    
    
    tic             % start the timer
    input = uint8(mean(double(imread(filenamelist((i+1),:))),3));       % read in the image which has to be correlated

    input_points_for(:,1)=reshape(input_points_x+displx(1,i),[],1);         % we reshape the input points to one row of values since this is the shape cpcorr will accept
    input_points_for(:,2)=reshape(input_points_y+disply(1,i),[],1);
    base_points_for(:,1)=reshape(base_points_x,[],1);
    base_points_for(:,2)=reshape(base_points_y,[],1);
    input_correl(:,:)=cpcorr(round(input_points_for), round(base_points_for), input, base);           % here we go and give all the markers and images to process to cpcorr.m which ic a function provided by the matlab image processing toolbox
    input_correl_x=input_correl(:,1);                                       % the results we get from cpcorr for the x-direction
    input_correl_y=input_correl(:,2);                                       % the results we get from cpcorr for the y-direction
    
    
    validx(:,i)=input_correl_x;                                                     % lets save the data
    savelinex=input_correl_x';
    dlmwrite('resultsimcorrx.txt', savelinex , 'delimiter', '\t', '-append');       % Here we save the result from each image; if you are desperately want to run this function with e.g. matlab 6.5 then you should comment this line out. If you do that the data will be saved at the end of the correlation step - good luck ;-)
    
    validy(:,i)=input_correl_y;
    saveliney=input_correl_y';
    dlmwrite('resultsimcorry.txt', saveliney , 'delimiter', '\t', '-append');
    
    waitbar(i/(r-1))                                                                        % update the waitbar
    
    % Update base and input points for cpcorr.m
    base_points_x=grid_x;
    base_points_y=grid_y;
    input_points_x=input_correl_x;
    input_points_y=input_correl_y;
    
    imshow(filenamelist(i+1,:))                     % update image
    hold on
    plot(grid_x,grid_y,'g+')                                % plot start position of raster
    plot(input_correl_x,input_correl_y,'r+')        % plot actual postition of raster
    hold off
    drawnow
    time(i)=toc;                                                 % take time
    estimatedtime=sum(time)/i*(r-1);            % estimate time to process
    title(['# Im.: ', num2str((r-1)),'; Proc. Im. #: ', num2str((i)),'; # Rasterp.:',num2str(row*col), '; Est. Time [s] ', num2str(round(estimatedtime)), ';  Elapsed Time [s] ', num2str(round(sum(time)))]);    % plot a title onto the image
    drawnow
    
end    

close(g)
close all

% save

save time.dat time -ascii -tabs
save validx.dat validx -ascii -tabs
save validy.dat validy -ascii -tabs
