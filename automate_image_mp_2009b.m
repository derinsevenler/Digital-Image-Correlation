function [validx,validy]=automate_image_mp_2009b(grid_x,grid_y,filenamelist,validx,validy);

%% Coarse grained multi core code to start actual image correlation
% Be aware that you will need the Parallel Computing Toolbox and at least
% Matlab 2009b, all earlier version do not work. Experience shows that
% logical processor cores as available within the Core i7 series of Intel
% do not accelerate the code much. this seems to be due to the highly
% optimized fft code of matlab and the extensive memory bandwidth needed.
%
% Last change by Chris
%
% Last revision: 9/8/10


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


[row,col]=size(base_points_x);      % this will determine the number of rasterpoints we have to run through
[r,c]=size(filenamelist);                   % this will determine the number of images we have to loop through

%% multicore part, here we start the matlabpool and the matlab workers
% after checking if there is a matlabpool available it starts a local one.
% If you want to limit the number of workers to the real core count (not
% the logical ones, you can add the number in line 94 (e.g. matlabpool 4)
clear ans
matlabpool size;
poolsize=ans;
if poolsize==0;
    matlabpool %start a new matlabpool
end
matlabpool size;
poolsize=ans;

grid_x=reshape(grid_x,[],1);
grid_y=reshape(grid_y,[],1);

numberperlab=length(grid_x)/poolsize;

% single program multiple data -  here we distribute the work between the
% available workers
spmd
if labindex<poolsize
    input_points_x=grid_x(1+(labindex-1)*fix(numberperlab):(labindex)*fix(numberperlab),1);
else
    input_points_x=grid_x(1+(labindex-1)*fix(numberperlab):length(grid_x),1);
end
    base_points_x=input_points_x;
end
spmd
if labindex<poolsize
    input_points_y=grid_y(1+(labindex-1)*fix(numberperlab):(labindex)*fix(numberperlab),1);
else
    input_points_y=grid_y(1+(labindex-1)*fix(numberperlab):length(grid_x),1);
end
    base_points_y=input_points_y;
end

%% Open new figure so previous ones (if open) are not overwritten
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

base = uint8(mean(double(imread(filenamelist(1,:))),3));            % read in the base image ( which is always  image number one. You might want to change that to improve correlation results in case the light conditions are changing during the experiment

for i=firstimage:(r-1)               % run through all images
    
    
    tic;             % start the timer
    input = uint8(mean(double(imread(filenamelist((i+1),:))),3));       % read in the image which has to be correlated
    
    % here we send of the image to the individual workers which will work
    % on their specific markers
    spmd
        input_points_for=[input_points_x input_points_y];
        base_points_for=[base_points_x base_points_y];
        % rounding the input and base points gives much better resolution -
        % seems to be a imcrop problem but I'm still investigating
        input_correl(:,:)=cpcorr(round(input_points_for),round(base_points_for), input, base);           % here we go and give all the markers and images to process to cpcorr.m which ic a function provided by the matlab image processing toolbox
        input_correl_x=input_correl(:,1);                                       % the results we get from cpcorr for the x-direction
        input_points_x=input_correl_x;
        input_correl_y=input_correl(:,2);                                       % the results we get from cpcorr for the y-direction
        input_points_y=input_correl_y;
    end
    
    % we need to collect all the data from our workers
    for j=1:numel(input_correl_x)
        lengthx(j,:)=size(input_correl_x{j});
    end
    bis=0;
    for j=1:numel(input_correl_x)
        von=bis+1;
        bis=bis+lengthx(j,1);
        validx_temp(von:bis,:)=input_correl_x{j};
        validy_temp(von:bis,:)=input_correl_y{j};
    end
    
    validx(:,i)=validx_temp;                                                     % lets save the data
    savelinex=validx_temp';
    dlmwrite('resultsimcorrx.txt', savelinex , 'delimiter', '\t', '-append');       % Here we save the result from each image; if you are desperately want to run this function with e.g. matlab 6.5 then you should comment this line out. If you do that the data will be saved at the end of the correlation step - good luck ;-)
    
    validy(:,i)=validy_temp;
    saveliney=validy_temp';
    dlmwrite('resultsimcorry.txt', saveliney , 'delimiter', '\t', '-append');
    
    waitbar(i/(r-1))   ;                                                                     % update the waitbar
    
    imshow(filenamelist(i+1,:))     ;                % update image
    hold on
    plot(grid_x,grid_y,'g+');                               % plot start position of raster
    plot(validx_temp,validy_temp,'r+');        % plot actual postition of raster
    hold off
    drawnow
    time(i)=toc;                                                 % take time
    estimatedtime=sum(time)/i*(r-1);           % estimate time to process
    title(['# Im.: ', num2str((r-1)),'; Proc. Im. #: ', num2str((i)),'; # Rasterp.:',num2str(row*col), '; Est. Time [s] ', num2str(estimatedtime), ';  Elapsed Time [s] ', num2str(sum(time))]);    % plot a title onto the image
    drawnow
    
end

close(g)
close all

% save
matlabpool close

save time.dat time -ascii -tabs
save validx.dat validx -ascii -tabs
save validy.dat validy -ascii -tabs
