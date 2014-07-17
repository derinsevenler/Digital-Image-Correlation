function [grid_x,grid_y]=grid_generator(FileNameBase,PathNameBase, grid_x, grid_y);

% Code to generate the DIC analysis grid
% Changed by Sven 
% Completely rewritten by Chris
% Programmed first by Dan and Rob 
% 
% Last revision: 05/09/10 Version 2.0

% The grid_generator function will help you create grids of markers.
% First you'll be asked for the base image that is used to define the grid
% which is typically your first image. Then you'll be asked if you want to 
% modify an old grid or create a new one. The dialog has different options
% allowing you to create a marker grid which is rectangular, circular, a
% line or contains only of two markers or delet markers from created grid.
% Every combination of them is also
% possible. You will be asked to click at the sites of interest and the
% markers will be plotted on top of your image. You can choose if you want
% to keep these markers or if you want to try again. If you keep them they
% will be saved and you'll come back to the main menu.
% It has to be noted that you can always generate your own marker 
% positions. Therefore the marker position in pixel has to be saved as a
% text based format where the x-position is saved as grid_x.dat and the
% y-position saved as grid_y.dat.
%

% Check if a grid is loaded if not new variables will be created
if exist('grid_x','var')~=1
    grid_x=[];
end

if exist('grid_y','var')~=1
    grid_y=[];
end


% Prompt user for base image if no image already assigned
if exist('FileNameBase','var')~=1
    [FileNameBase,PathNameBase,FilterIndex] = uigetfile( ...
    {'*.bmp;*.tif;*.jpg;*.TIF;*.BMP;*.JPG','Image files (*.bmp,*.tif,*.jpg)';'*.*',  'All Files (*.*)'}, ...
    'Open base image for grid creation');
else 
    FilterIndex=1;
end
FilterIndex
FileNameBase
PathNameBase
% Check if an image is chosen, if true go to directory an load image, call
% gridtypeselection, if false end program
if FilterIndex~=0 
     cd(PathNameBase)   
     im_grid = imread(FileNameBase); 
     [grid_x,grid_y,FileNameBase,PathNameBase] = gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,1); %call gridselection
else %if not the empty variables will be cleared
     clear  FileNameBase PathNameBase FilterIndex
	 disp('No Image is chosen')
end

close all
end
%% Main Menu
function [grid_x,grid_y,FileNameBase,PathNameBase] = gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,loadgrid)

% Code to select which kind of grid should be added and to display the actual grid
% Decide which type of grid you want to create
% In this area you should select the type of grid you want to add or create
%
hold off
imshow(im_grid,'InitialMagnification',100); %show chosen Image
%------------------------
% Load old grid?
% 
if loadgrid
loadoldgrid=menu(sprintf('Load old grid?'),...
'Yes','No');

if loadoldgrid==1
	[gridxname,Pathgridx] = uigetfile('*.dat','Open grid_x.dat'); %load grid_x
		if gridxname==0
			disp('You did not select a file!')
        end
    cd(Pathgridx);
    grid_x=importdata(gridxname,'\t');
	
	[gridyname,Pathgridy] = uigetfile('*.dat','Open grid_y.dat');%load grid_y
		if gridyname==0
			disp('You did not select a file!')
        end
    cd(Pathgridy);
    grid_y=importdata(gridyname,'\t');
    
    
end
end
hold on %plot old grid
plot(grid_x, grid_y,'+r')
hold off
%------------------------
% Grid selection
%
gridselection = menu(sprintf('Gridgenerator Main Menu'),...
    'Rectangular','Circular','Two Markers','Line','Remove Markers','END');


    % Rectangular
    if gridselection==1
        [grid_x,grid_y,FileNameBase,PathNameBase] = rect_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
        
    end
    % Circular
    if gridselection==2
        [grid_x,grid_y,FileNameBase,PathNameBase] = circ_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
        
    end
    % Two Markers
    if gridselection==3
        [grid_x,grid_y,FileNameBase,PathNameBase] = twop_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
        
    end
    % Line
    if gridselection==4
        [grid_x,grid_y,FileNameBase,PathNameBase] = line_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
        
    end
    % Remove Markers
    if gridselection==5
    [grid_x,grid_y,FileNameBase,PathNameBase] = removepoints(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
    
    end
    % END
    if gridselection==6
    end
end
%% Select a rect area
function [grid_x,grid_y,FileNameBase,PathNameBase] = rect_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y)

% Function to select a rectangular grid and to add these to an existing one
% wirtten by Chris
%

title(sprintf('Define the region of interest.  Pick (single click) a point in the LOWER LEFT region of the gage section.\n  Do the same for a point in the UPPER RIGHT portion of the gage section.'))

[x(1,1),y(1,1)]=ginput(1);
hold on
plot(x(1,1),y(1,1),'+b')

[x(2,1),y(2,1)]=ginput(1);
hold on
plot(x(2,1),y(2,1),'+b')

drawnow

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

lowerline=[xmin ymin; xmax ymin];
upperline=[xmin ymax; xmax ymax];
leftline=[xmin ymin; xmin ymax];
rightline=[xmax ymin; xmax ymax];

plot(lowerline(:,1),lowerline(:,2),'-b')
plot(upperline(:,1),upperline(:,2),'-b')
plot(leftline(:,1),leftline(:,2),'-b')
plot(rightline(:,1),rightline(:,2),'-b')

% closereq

cd(PathNameBase)

% Prompt user for grid spacing/resolution
prompt = {'Enter horizontal (x) resolution for image analysis [pixels]:', ...
        'Enter vertical (y) resolution for image analysis [pixels]:'};
dlg_title = 'Input for grid creation';
num_lines= 1;
def     = {'15','15'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
xspacing = str2double(cell2mat(answer(1,1)));
yspacing = str2double(cell2mat(answer(2,1)));

% Round xmin,xmax and ymin,ymax "up" based on selected spacing
numXelem = ceil((xmax-xmin)/xspacing)-1;
numYelem = ceil((ymax-ymin)/yspacing)-1;

xmin_new = (xmax+xmin)/2-((numXelem/2)*xspacing);
xmax_new = (xmax+xmin)/2+((numXelem/2)*xspacing);
ymin_new = (ymax+ymin)/2-((numYelem/2)*yspacing);
ymax_new = (ymax+ymin)/2+((numYelem/2)*yspacing);

% Create the analysis grid and show user
[x,y] = meshgrid(xmin_new:xspacing:xmax_new,ymin_new:yspacing:ymax_new);
[rows columns] = size(x);
%zdummy = 200.*ones(rows,columns);
imshow(FileNameBase)
title(['Selected grid has ',num2str(rows*columns), ' rasterpoints'])    % plot a title onto the image
hold on;
plot(grid_x,grid_y,'+r')
plot(x,y,'+b')


% Do you want to keep/add the grid?
confirmselection = menu(sprintf('Do you want to use this grid?'),...
    'Yes','No, try again','Go back to Main Menu');

    % Yes
    if confirmselection==1
        % Save settings and grid files in the image directory for visualization/plotting later
        x=reshape(x,[],1);
        y=reshape(y,[],1);
        grid_x=[grid_x;x];
        grid_y=[grid_y;y];
        save settings.dat xspacing yspacing xmin_new xmax_new ymin_new ymax_new -ascii -tabs
        save grid_x.dat grid_x -ascii -tabs
        save grid_y.dat grid_y -ascii -tabs

        close all
        hold off
        gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
    end

    % No, try again
    if confirmselection==2
        close all
        hold off
        imshow(im_grid,'truesize');
        hold on;
        plot(grid_x,grid_y,'+r');
        hold off;
        rect_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
    end
    
    % Go back to Main Menu
    if confirmselection==3
        close all
        hold off
        gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
    end
end
%% Select a circular area
function [grid_x,grid_y,FileNameBase,PathNameBase] = circ_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y)

title(sprintf('Pick three points on the circle in clockwise order at the upper boundary of the sample.') )

[x(1,1),y(1,1)]=ginput(1);
hold on
plot(x(1,1),y(1,1),'+g')

[x(2,1),y(2,1)]=ginput(1);
plot(x(2,1),y(2,1),'+g')

[x(3,1),y(3,1)]=ginput(1);
plot(x(3,1),y(3,1),'+g')

xnew=x;
ynew=y;

% Calculate center between the 3 sorted points and the normal slope of the vectors
slope12=-1/((ynew(2,1)-ynew(1,1))/(xnew(2,1)-xnew(1,1)));
slope23=-1/((ynew(3,1)-ynew(2,1))/(xnew(3,1)-xnew(2,1)));
center12(1,1)=(xnew(2,1)-xnew(1,1))/2+xnew(1,1);
center12(1,2)=(ynew(2,1)-ynew(1,1))/2+ynew(1,1);
center23(1,1)=(xnew(3,1)-xnew(2,1))/2+xnew(2,1);
center23(1,2)=(ynew(3,1)-ynew(2,1))/2+ynew(2,1);



% Calculate the crossing point of the two vectors
achsenabschnitt1=center12(1,2)-center12(1,1)*slope12;
achsenabschnitt2=center23(1,2)-center23(1,1)*slope23;
xcross=(achsenabschnitt2-achsenabschnitt1)/(slope12-slope23);
ycross=slope12*xcross+achsenabschnitt1;
plot(xcross,ycross,'or')

% Calculate radius 
R=sqrt((xcross-xnew(1,1))*(xcross-xnew(1,1))+(ycross-ynew(1,1))*(ycross-ynew(1,1)));

% Calculate angle between vectors
xvector=[1;0];
x1vec(1,1)=xnew(1,1)-xcross;x1vec(2,1)=ynew(1,1)-ycross;
x3vec(1,1)=xnew(3,1)-xcross;x3vec(2,1)=ynew(3,1)-ycross;
alpha13=acos((dot(x1vec,x3vec))/(sqrt(x1vec'*x1vec)*sqrt(x3vec'*x3vec)))*180/pi;
alpha03=acos((dot(xvector,x3vec))/(sqrt(xvector'*xvector)*sqrt(x3vec'*x3vec)))*180/pi;
totalangle=alpha13;
maxangle=alpha03;
angldiv=abs(round(totalangle))*10;
anglstep=(totalangle/angldiv);
anglall(1:angldiv+1)=maxangle+anglstep*(1:angldiv+1)-anglstep;
xcircle(1:angldiv+1)=xcross+R*cos(-anglall(1:angldiv+1)/180*pi);
ycircle(1:angldiv+1)=ycross+R*sin(-anglall(1:angldiv+1)/180*pi);
plot(xcircle,ycircle,'-b')
drawnow

title(['Segment of circle spreads over ',num2str(totalangle),'°'])


% Accept the chosen circle, try again or give up 

confirmcircselection = menu(sprintf('Do you want to use this circle as basis?'),...
    'Yes','No, try again','Go back to grid-type selection');
    
    % No, try again
    if confirmcircselection==2
        close all
        imshow(im_grid,'truesize');
        hold on
        plot(grid_x,grid_y,'+r');
        hold off
        circ_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
    end
    
    % Go back to grid-type selection
    if confirmcircselection==3
        close all
        gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
    end

    % Yes
    if confirmcircselection==1

        prompt = {'Enter the number of intersections between markers on the circle:'};
        dlg_title = 'Input for grid creation';
        num_lines= 1;
        def     = {'30'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        angldiv = str2double(cell2mat(answer(1,1)));

        anglstep=(totalangle/angldiv);
        anglall(1:angldiv+1)=maxangle+anglstep*(1:angldiv+1)-anglstep;

        markerxpos(1:angldiv+1)=xcross+R*cos(-anglall(1:angldiv+1)/180*pi);
        markerypos(1:angldiv+1)=ycross+R*sin(-anglall(1:angldiv+1)/180*pi);

        plot(markerxpos,markerypos,'ob');

        % Pick the lower bound in the image
        title(sprintf('Pick three points lying on the circle in clockwise order. The first and last one define the width of the raster') )

        [x(4,1),y(4,1)]=ginput(1);
        hold on
        plot(x(1,1),y(1,1),'+r')

        lowboundx=x(4,1);
        lowboundy=y(4,1);

        R2=sqrt((xcross-lowboundx(1,1))*(xcross-lowboundx(1,1))+(ycross-lowboundy(1,1))*(ycross-lowboundy(1,1)));
        markerxposlb(1:angldiv+1)=xcross+R2*cos(-anglall(1:angldiv+1)/180*pi);
        markeryposlb(1:angldiv+1)=ycross+R2*sin(-anglall(1:angldiv+1)/180*pi);

        plot(markerxposlb,markeryposlb,'ob');

        prompt = {'Enter the number of intersections between the upper and lower bound:'};
        dlg_title = 'Input for grid creation';
        num_lines= 1;
        def     = {'5'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        Rdiv = str2double(cell2mat(answer(1,1)));

        Rstep=((R-R2)/Rdiv);
        Rall(1:Rdiv+1)=R2+Rstep*(1:Rdiv+1)-Rstep;

        x=ones(Rdiv+1,angldiv+1)*xcross;
        y=ones(Rdiv+1,angldiv+1)*ycross;
        x=x+Rall'*cos(-anglall(1:angldiv+1)/180*pi);
        y=y+Rall'*sin(-anglall(1:angldiv+1)/180*pi);

        close all
        imshow(im_grid,'truesize');
        hold on
        plot(grid_x,grid_y,'+r')    
        plot(x,y,'.b')    

        title(['Selected grid has ',num2str(angldiv*Rdiv), ' rasterpoints'])    % plot a title onto the image


         % Do you want to keep/add the grid?
        confirmselection = menu(sprintf('Do you want to use this grid?'),...
        'Yes','No, try again','Go back to Main Menu');

        % Yes
        if confirmselection==1
            % Save settings and grid files in the image directory for visualization/plotting later
            x=reshape(x,[],1);
            y=reshape(y,[],1);
            grid_x=[grid_x;x];
            grid_y=[grid_y;y];
            save grid_x.dat grid_x -ascii -tabs
            save grid_y.dat grid_y -ascii -tabs
            close all
            hold off
            gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
        end

        % No, try again
        if confirmselection==2
            close all
            hold off
            imshow(im_grid,'truesize');
            hold on;
            plot(grid_x,grid_y,'+r');
            hold off;
            circ_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
        end

        % Go back to Main Menu
        if confirmselection==3
            close all
            hold off
            gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
        end
    end
end
%% Select 2 Points
function [grid_x,grid_y,FileNameBase,PathNameBase] = twop_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y)

title(sprintf('Pick two points on the sample.') )

[x(1,1),y(1,1)]=ginput(1);
hold on
plot(x(1,1),y(1,1),'+b')

[x(2,1),y(2,1)]=ginput(1);
plot(x(2,1),y(2,1),'+b')

% Do you want to keep/add the grid?
confirmselection = menu(sprintf('Do you want to use this grid?'),...
    'Yes','No, try again','Go back to Main Menu');

    % Yes
    if confirmselection==1
        % Save settings and grid files in the image directory for visualization/plotting later
        x=reshape(x,[],1);
        y=reshape(y,[],1);
        grid_x=[grid_x;x];
        grid_y=[grid_y;y];
        save grid_x.dat grid_x -ascii -tabs
        save grid_y.dat grid_y -ascii -tabs
        close all
        hold off
        gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
    end

    % No, try again
    if confirmselection==2
        close all
        hold off
        imshow(im_grid,'truesize');
        hold on;
        plot(grid_x,grid_y,'+r');
        hold off;
        twop_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
    end
    
    % Go back to Main Menu
    if confirmselection==3
        close all
        hold off
        gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
    end
end
%% Select a line
function [grid_x,grid_y,FileNameBase,PathNameBase] = line_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y)

title(sprintf('Pick two points on the sample.') )

[x(1,1),y(1,1)]=ginput(1);
hold on
plot(x(1,1),y(1,1),'+b')

[x(2,1),y(2,1)]=ginput(1);
plot(x(2,1),y(2,1),'+b')

lineslope=(y(2,1)-y(1,1))/(x(2,1)-x(1,1));
intersecty=y(1,1)-lineslope*x(1,1);
ycalc=zeros(2,1);
ycalc=lineslope*x+intersecty;
plot(x(:,1),ycalc(:,1),'-b')


prompt = {'Enter the number of intersections between markers on the line:'};
dlg_title = 'Input for grid creation';
num_lines= 1;
def     = {'30'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
linediv = str2num(cell2mat(answer(1,1)));
linestep=((max(x)-min(x))/linediv);
x(1:linediv+1)=min(x)+linestep*(1:linediv+1)-linestep;
y=lineslope*x+intersecty;

plot(x,y,'ob')
title(['Selected grid has ',num2str(linediv), ' rasterpoints'])    % plot a title onto the image


% Do you want to keep/add the grid?
confirmselection = menu(sprintf('Do you want to use this grid?'),...
    'Yes','No, try again','Go back to Main Menu');

    % Yes
    if confirmselection==1
        % Save settings and grid files in the image directory for visualization/plotting later
        x=reshape(x,[],1);
        y=reshape(y,[],1);
        grid_x=[grid_x;x];
        grid_y=[grid_y;y];
        save grid_x.dat grid_x -ascii -tabs
        save grid_y.dat grid_y -ascii -tabs
        close all
        hold off
        gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
    end

    % No, try again
    if confirmselection==2
        close all
        hold off
        imshow(im_grid,'truesize');
        hold on;
        plot(grid_x,grid_y,'+r');
        hold off;
        line_grid(FileNameBase, PathNameBase, im_grid, grid_x, grid_y);
    end
    
    % Go back to Main Menu
    if confirmselection==3
        close all
        hold off
        gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
    end
end
%% Select Points to remove
function [grid_x,grid_y,FileNameBase,PathNameBase] = removepoints(FileNameBase, PathNameBase, im_grid, grid_x, grid_y)

% Delete some markers
%

%create working copy of the grid

grid_xtemp=grid_x;
grid_ytemp=grid_y;

imshow(im_grid,'InitialMagnification',100); %show chosen Image
hold on
plot(grid_x, grid_y,'ob');
hold off,

title(sprintf('Define the region of interest.  \n  All points inside that region will be deleted'))

        [xdel,ydel]=ginput(2);
        x(1,1) = xdel(1);
        x(1,2) = xdel(2);
        y(1,1) = ydel(2);
        y(1,2) = ydel(1);

        deletepoints=find(grid_x>min(x) & grid_x<max(x) & grid_y<max(y) & grid_y>min(y));
        
        grid_xtemp(deletepoints,:)=[];
        grid_ytemp(deletepoints,:)=[];

        imshow(im_grid,'InitialMagnification',100); %show chosen Image
        hold on
        plot(grid_xtemp, grid_ytemp,'ob');
        hold off,

        
        % delete point permanently?
        keepchanges = menu(sprintf('Do you want to delete these markers permanently?'),'Yes','No');
        if keepchanges==1
            grid_x=grid_xtemp;
            grid_y=grid_ytemp;
            save grid_x.dat grid_x -ascii -tabs
            save grid_y.dat grid_y -ascii -tabs
            close all
            hold off
            gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
        end
        
        if keepchanges==2
            gridtypeselection(FileNameBase, PathNameBase, im_grid, grid_x, grid_y,0);
        end
        
        
end