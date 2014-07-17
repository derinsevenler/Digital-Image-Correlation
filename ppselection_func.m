% written by Chris

function [rasterx, rastery, validx,validy,x,y]=ppselection_func(validx,validy,x,y);

% Code to analyze the displacement data (contained by validx and validy)
% Programmed by Chris
% Last revision: 8/24/06

% Prompt user for displacement data
if exist('validx')==0
    [validx,Pathvalidx] = uigetfile('*.mat','Open validx.mat');
    if validx==0
        return;
    end
    cd(Pathvalidx);
    validx=importdata(validx,'\t');
end
if exist('validy')==0
    [validy,Pathvalidy] = uigetfile('*.mat','Open validy.mat');
    if validy==0
        return;
    end
    cd(Pathvalidy);
    validy=importdata(validy,'\t');
end

% Checking for plot orientation and give standard orientations
if exist('x')==0
    x=1;
    y=2;
end
if x~1|2|3
    x=1;
    y=2;
end

validxbackup=validx;
validybackup=validy;

% Choose an image
[looppoints loopimages]=size(validx);
selectedimage=0;
prompt = {'From which image do you want to select the view?'};
dlg_title = 'Marker selection';
num_lines= 1;
if selectedimage==0
    defaultimage=loopimages;
end
if selectedimage~0
    defaultimage=selectedimage;
end
def     = {num2str(defaultimage)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
selectedimage = str2num(cell2mat(answer(1,1)));
if selectedimage>loopimages
    selectedimage=loopimages;
end
if selectedimage<1
    selectedimage=1;
end
rasterx=0;
rastery=0;
% Call the selection tool
[rasterx, rastery, validx,validy,x,y]=gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);

return

%-------------------------------
%
% Decide which type of raster you want to analyze

function [rasterx, rastery, validx,validy,x,y] = gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);

gridselection = menu(sprintf('Which type of grid do you want to use'),...
    'Two Markers','Rectangular','Two Rectangles of Markers','Circular','Line','Change view','Cancel');

if gridselection==1
    [validx,validy,rasterx, rastery,selectedimage,x,y]=twop_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    exist('validx');
    exist('validy');
    exist('rasterx');
    exist('rastery');
    exist('x');
    exist('y');
    return
end

if gridselection==2
    [validx,validy,rasterx, rastery,selectedimage,x,y]=rect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    exist('validx');
    exist('validy');
    exist('rasterx');
    exist('rastery');
    exist('x');
    exist('y');
    return
end

if gridselection==5
    [validx,validy,rasterx, rastery,selectedimage,x,y]=line_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    exist('validx');
    exist('validy');
    exist('rasterx');
    exist('rastery');
    exist('x');
    exist('y');
    return
end

if gridselection==3
    [validx,validy,rasterx, rastery,selectedimage,x,y]=tworect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    exist('validx');
    exist('validy');
    exist('rasterx');
    exist('rastery');
    exist('x');
    exist('y');
    return
end

if gridselection==4
    [validx,validy,rasterx, rastery,selectedimage,x,y]=circ_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    exist('validx');
    exist('validy');
    exist('rasterx');
    exist('rastery');
    exist('x');
    exist('y');
    return
end

if gridselection==6
    [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    exist('validx');
    exist('validy');
    exist('rasterx');
    exist('rastery');
    exist('x');
    exist('y');
end

if gridselection==7
    close all
    exist('validx');
    exist('validy');
    exist('rasterx');
    exist('rastery');
    exist('x');
    exist('y');
    if exist('validx')==0
        validx=0;
        validy=0;
    end
    if exist('rasterx')==0
        rasterx=0;
        rastery=0;
    end
    if exist('x')==0
        x=0;
        y=0;
    end
end

if exist('validx')==0
    validx=0
    validy=0
end
if exist('rasterx')==0
    rasterx=0
    rastery=0
end
if exist('x')==0
    x=0
    y=0
end

return

%-------------------------------
%
% Change the view

function [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);

[looppoints loopimages]=size(validx);
if x==1
    pos=validx;
end
if x==2
    pos=validy;
end
if x==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    pos=validx-validxfirst;
end
if x==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    pos=validy-validyfirst;
end
if y==1
    displ=validx;
end
if y==2
    displ=validy;
end
if y==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    displ=validx-validxfirst;
end
if y==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    displ=validy-validyfirst;
end

h=figure;
plot(pos(:,selectedimage),displ(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','g')

title(sprintf('View.') )

changeviewselection = menu(sprintf('Do you want to change the coordinate system to select markers?'),...
    'x-position vs. y-position','x-position vs. x-displacement','y-position vs. x-displacement',...
    'x-position vs. y-displacement','y-position vs. y-displacement','Change image #','Go back to grid-type selection');

if changeviewselection==1
    close all
    x=1;
    y=2;
    [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end
if changeviewselection==2
    close all
    x=1;
    y=3;
    [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end
if changeviewselection==3
    close all
    x=2;
    y=3;
    [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end
if changeviewselection==4
    close all
    x=1;
    y=4;
    [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end
if changeviewselection==5
    close all
    x=2;
    y=4;
    [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end
if changeviewselection==6
    prompt = {'From which image do you want to select the view?'};
    dlg_title = 'Marker selection';
    num_lines= 1;
    if selectedimage==0
        defaultimage=loopimages;
    end
    if selectedimage~0
        defaultimage=selectedimage;
    end
    def     = {num2str(defaultimage)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    selectedimage = str2num(cell2mat(answer(1,1)));
    if selectedimage>loopimages
        selectedimage=loopimages;
    end
    if selectedimage<1
        selectedimage=1;
    end
    [validx,validy,x,y,selectedimage,rasterx,rastery] = change_view_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end
if changeviewselection==7
    close all
    [rasterx, rastery, validx,validy,x,y]=gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);
end

%-------------------------------
%
% Define two rectangles and add them to one marker array

function [validx,validy,rasterx,rastery,selectedimage,x,y] = tworect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);

[validx,validy,rasterx1, rastery1,selectedimage]=rect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
[validx,validy,rasterx2, rastery2,selectedimage]=rect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);

% if size(rasterx1)~size(rasterx2)
%     close all
%     rasterx=0;rastery=0;
%     return
% end

rasterx=[rasterx1; rasterx2];
rastery=[rastery1; rastery2];

if x==1
    pos1=validx;pos2=rasterx;
end
if x==2
    pos1=validy;pos2=rasterx;
end
if x==3
    validxfirst=zeros(size(validx));rasterxfirst=zeros(size(rasterx));
    validxfirst=validx(:,1)*ones(1,loopimages);rasterxfirst=rasterx(:,1)*ones(1,loopimages);
    pos1=validx-validxfirst;pos2=rasterx-rasterxfirst;
end
if x==4
    validyfirst=zeros(size(validy));rasteryfirst=zeros(size(rastery));
    validyfirst=validy(:,1)*ones(1,loopimages);rasteryfirst=rastery(:,1)*ones(1,loopimages);
    pos1=validy-validyfirst;pos2=rastery-rasteryfirst;
end
if y==1
    displ1=validx;displ2=rasterx;
end
if y==2
    displ1=validy;displ2=rastery;
end
if y==3
    validyfirst=zeros(size(validy));rasteryfirst=zeros(size(rastery));
    validyfirst=validy(:,1)*ones(1,loopimages);rasteryfirst=rastery(:,1)*ones(1,loopimages);
    displ1=validy-validyfirst;displ2=rastery-rasteryfirst;
end

[looppoints loopimages]=size(validx);
defaultimage=loopimages;

prompt = {'From which image do you want to select the rectangles?'};
dlg_title = 'Marker selection';
num_lines= 1;
if selectedimage==0
    defaultimage=loopimages;
end
if selectedimage~0
    defaultimage=selectedimage;
end
def     = {num2str(defaultimage)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
selectedimage = str2num(cell2mat(answer(1,1)));
if selectedimage>loopimages
    selectedimage=loopimages;
end
if selectedimage<1
    selectedimage=1;
end

h=figure;
plot(pos1(:,selectedimage),displ1(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
hold on
plot(pos2(:,selectedimage),displ2(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','r')

% Accept the chosen markers, try again or give up 

confirmcircselection = menu(sprintf('Do you want to use these markers?'),...
    'Yes','No, try again','Go back to grid-type selection');

if confirmcircselection==2
    close all
    [validx,validy,rasterx,rastery,selectedimage,x,y] = tworect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end

if confirmcircselection==3
    close all
    [rasterx, rastery,validx,validy, x,y] = gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);
end

if confirmcircselection==1
    close all
end

%-------------------------------
%
% Define line and find markers

function [validx,validy,rasterx,rastery,selectedimage,x,y]=line_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);

[looppoints loopimages]=size(validx);
defaultimage=loopimages;

if x==1
    pos=validx;
end
if x==2
    pos=validy;
end
if x==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    pos=validx-validxfirst;
end
if x==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    pos=validy-validyfirst;
end
if y==1
    displ=validx;
end
if y==2
    displ=validy;
end
if y==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    displ=validx-validxfirst;
end
if y==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    displ=validy-validyfirst;
end

prompt = {'From which image do you want to select markers?'};
dlg_title = 'Marker selection';
num_lines= 1;
if selectedimage==0
    defaultimage=loopimages;
end
if selectedimage~0
    defaultimage=selectedimage;
end
def     = {num2str(defaultimage)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
selectedimage = str2num(cell2mat(answer(1,1)));
if selectedimage>loopimages
    selectedimage=loopimages;
end
if selectedimage<1
    selectedimage=1;
end

h=figure;
plot(pos(:,selectedimage),displ(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','g')

title(sprintf('Pick two points on the sample.') )

[xpos(1,1),ypos(1,1)]=ginput(1);
hold on
plot(xpos(1,1),ypos(1,1),'+g')

[xpos(2,1),ypos(2,1)]=ginput(1);
plot(xpos(2,1),ypos(2,1),'+g')

centerpoint=[xpos(1,1)+(xpos(2,1)-xpos(1,1))/2; ypos(1,1)+(ypos(2,1)-ypos(1,1))/2];
plot(centerpoint(1,1),centerpoint(2,1),'+b')

linelength=sqrt((xpos(2,1)-xpos(1,1))*(xpos(2,1)-xpos(1,1))+(ypos(2,1)-ypos(1,1))*(ypos(2,1)-ypos(1,1)));
lineslope=(ypos(2,1)-ypos(1,1))/(xpos(2,1)-xpos(1,1));
intersecty=ypos(1,1)-lineslope*xpos(1,1);
ycalc=zeros(2,1);
ycalc=lineslope*xpos+intersecty;
plot(xpos(:,1),ycalc(:,1),'-b')
intercept=[0; centerpoint(2,1)-centerpoint(1,1)*(-1/lineslope)];

distancefromline=(abs((xpos(2,1)-xpos(1,1))*(ypos(1,1)-displ(:,selectedimage))-(xpos(1,1)-pos(:,selectedimage))*(ypos(2,1)-ypos(1,1))))/sqrt((xpos(2,1)-xpos(1,1))*(xpos(2,1)-xpos(1,1))+(ypos(2,1)-ypos(1,1))*(ypos(2,1)-ypos(1,1)));
distancefromcenterpoint=(abs((intercept(1,1)-centerpoint(1,1))*(centerpoint(2,1)-displ(:,selectedimage))-(centerpoint(1,1)-pos(:,selectedimage))*(intercept(2,1)-centerpoint(2,1))))/sqrt((intercept(1,1)-centerpoint(1,1))*(intercept(1,1)-centerpoint(1,1))+(intercept(2,1)-centerpoint(2,1))*(intercept(2,1)-centerpoint(2,1)));

linewidthquestion=0;
linewidth=20;

while linewidthquestion==0
    prompt = {'Enter the width of the line in [pixel]:'};
    dlg_title = 'Input for grid creation';
    num_lines= 1;
    def     = {num2str(linewidth)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    linewidth = str2num(cell2mat(answer(1,1)))
    
    selectpoints=find(distancefromline<linewidth & distancefromcenterpoint<linelength/2);
    plot(pos(selectpoints,selectedimage),displ(selectpoints,selectedimage),'.r')
    drawnow
    
    confirmlineselection = menu(sprintf('Do you want to use these markers?'),...
        'Yes','No, try different linewidth','No, different line','Go back to grid-type selection');
    
    if confirmlineselection==1
        close all
        linewidthquestion=1;
        rasterx=[validx(selectpoints,:);validx(selectpoints,:)];
        rastery=[validy(selectpoints,:);validy(selectpoints,:)];
        close all
    end
    if confirmlineselection==2
        hold off
        plot(pos(:,selectedimage),displ(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
        hold on
        plot(xpos(:,1),ycalc(:,1),'-b')
        plot(xpos(1,1),ypos(1,1),'+g')
        plot(xpos(2,1),ypos(2,1),'+g')
    end
    if confirmlineselection==3
        close all
        linewidthquestion=1;
        [validx,validy,rasterx, rastery,selectedimage,x,y]=line_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    end
    
    if confirmlineselection==4
        close all
        linewidthquestion=1;
        [rasterx, rastery,validx,validy, x,y] = gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);
    end
    
end

%-------------------------------
%
% Choose two markers

function [validx,validy,rasterx, rastery,selectedimage,x,y]=twop_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);

[looppoints loopimages]=size(validx);
defaultimage=loopimages;

if x==1
    pos=validx;
end
if x==2
    pos=validy;
end
if x==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    pos=validx-validxfirst;
end
if x==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    pos=validy-validyfirst;
end
if y==1
    displ=validx;
end
if y==2
    displ=validy;
end
if y==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    displ=validx-validxfirst;
end
if y==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    displ=validy-validyfirst;
end
prompt = {'From which image do you want to select markers?'};
dlg_title = 'Marker selection';
num_lines= 1;
if selectedimage==0
    defaultimage=loopimages;
end
if selectedimage~0
    defaultimage=selectedimage;
end
def     = {num2str(defaultimage)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
selectedimage = str2num(cell2mat(answer(1,1)));
if selectedimage>loopimages
    selectedimage=loopimages;
end
if selectedimage<1
    selectedimage=1;
end

h=figure;
plot(pos(:,selectedimage),displ(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','g')

title(sprintf('Pick two markers.') )

[xpos(1,1),ypos(1,1)]=ginput(1);
whereisthispoint=abs(pos(:,selectedimage)-xpos(1,1))+abs(displ(:,selectedimage)-ypos(1,1));
selectedpoint1=find(whereisthispoint==min(whereisthispoint));
hold on
plot(pos(selectedpoint1,selectedimage),displ(selectedpoint1,selectedimage),'+r')

[xpos(2,1),ypos(2,1)]=ginput(1);
whereisthispoint=abs(pos(:,selectedimage)-xpos(2,1))+abs(displ(:,selectedimage)-ypos(2,1));
selectedpoint2=find(whereisthispoint==min(whereisthispoint));
plot(pos(selectedpoint2,selectedimage),displ(selectedpoint2,selectedimage),'+r')


% Accept the chosen markers, try again or give up 

confirmcircselection = menu(sprintf('Do you want to use these two markers?'),...
    'Yes','No, try again','Go back to grid-type selection');

if confirmcircselection==2
    close all
    hold off
    [validx,validy,rasterx,rastery,selectedimage,x,y]=twop_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end

if confirmcircselection==3
    close all
    [rasterx,rastery,validx,validy,x,y]=gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);
end

if confirmcircselection==1
    rasterx=[validx(selectedpoint1,:);validx(selectedpoint2,:)];
    rastery=[validy(selectedpoint1,:);validy(selectedpoint2,:)];
    close all
end

%-------------------------------
%
%Circular selection

function [validx,validy,rasterx,rastery,selectedimage,x,y] = circ_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);

if x==1
    pos=validx;
end
if x==2
    pos=validy;
end
if x==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    pos=validx-validxfirst;
end
if x==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    pos=validy-validyfirst;
end
if y==1
    displ=validx;
end
if y==2
    displ=validy;
end
if y==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    displ=validx-validxfirst;
end
if y==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    displ=validy-validyfirst;
end
h=figure;
plot(pos(:,selectedimage),displ(:,selectedimage),'.b')
% axis equal

title(sprintf('Pick three points on the circle in clockwise order with the highest radius.') )

[xpos(1,1),ypos(1,1)]=ginput(1);
hold on
plot(xpos(1,1),ypos(1,1),'+g')

[xpos(2,1),ypos(2,1)]=ginput(1);
plot(xpos(2,1),ypos(2,1),'+g')

[xpos(3,1),ypos(3,1)]=ginput(1);
plot(xpos(3,1),ypos(3,1),'+g')

% Calculate center between the 3 sorted points and the normal slope of the vectors
slope12=-1/((ypos(2,1)-ypos(1,1))/(xpos(2,1)-xpos(1,1)));
slope23=-1/((ypos(3,1)-ypos(2,1))/(xpos(3,1)-xpos(2,1)));
center12(1,1)=(xpos(2,1)-xpos(1,1))/2+xpos(1,1);
center12(1,2)=(ypos(2,1)-ypos(1,1))/2+ypos(1,1);
center23(1,1)=(xpos(3,1)-xpos(2,1))/2+xpos(2,1);
center23(1,2)=(ypos(3,1)-ypos(2,1))/2+ypos(2,1);
plot(center12(1,1),center12(1,2),'+b')
plot(center23(1,1),center23(1,2),'+b')

if slope12==slope23
    [validx,validy,rasterx,rastery,selectedimage,x,y] = circ_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery)
end

% Calculate the crossing point of the two vectors
achsenabschnitt1=center12(1,2)-center12(1,1)*slope12;
achsenabschnitt2=center23(1,2)-center23(1,1)*slope23;
xcross=(achsenabschnitt2-achsenabschnitt1)/(slope12-slope23);
ycross=slope12*xcross+achsenabschnitt1;
xdata=min(xpos):xcross;
ydata1=achsenabschnitt1+slope12*xdata;
ydata2=achsenabschnitt2+slope23*xdata;

% Calculate radius
R=sqrt((xcross-xpos(1,1))*(xcross-xpos(1,1))+(ycross-ypos(1,1))*(ycross-ypos(1,1)))

% Calculate angle between vectors
xvector=[1;0];
x1vec(1,1)=xpos(1,1)-xcross;x1vec(2,1)=ypos(1,1)-ycross;
x3vec(1,1)=xpos(3,1)-xcross;x3vec(2,1)=ypos(3,1)-ycross;
alpha13=acos((dot(x1vec,x3vec))/(sqrt(x1vec'*x1vec)*sqrt(x3vec'*x3vec)))*180/pi;
alpha01=acos((dot(xvector,x1vec))/(sqrt(x1vec'*x1vec)*sqrt(xvector'*xvector)))*180/pi;
if ypos(1,1)<ycross
    alpha01=alpha01*(-1)+360;
end
alpha03=acos((dot(xvector,x3vec))/(sqrt(xvector'*xvector)*sqrt(x3vec'*x3vec)))*180/pi;
if ypos(3,1)<ycross
    alpha03=alpha03*(-1)+360;
end
totalangle=alpha13
minangle=alpha01
maxangle=alpha03

angldiv=abs(round(totalangle))*10;
anglstep=(totalangle/angldiv);
anglall(1:angldiv+1)=minangle-anglstep*(1:angldiv+1)-anglstep;
xcircle(1:angldiv+1)=xcross+R*cos(anglall(1:angldiv+1)/180*pi);
ycircle(1:angldiv+1)=ycross+R*sin(anglall(1:angldiv+1)/180*pi);
plot(xcircle,ycircle,'-r');
drawnow

% Accept the chosen circle, try again or give up 
confirmcircselection = menu(sprintf('Do you want to use this circle as basis?'),...
    'Yes','No, try again','Go back to grid-type selection');
if confirmcircselection==2
    close all
    hold off
    [validx,validy,rasterx,rastery,selectedimage,x,y] = circ_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end

if confirmcircselection==3
    close all
    [rasterx,rastery,validx,validy,x,y]=gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);
end

if confirmcircselection==1
    
    % Pick the lower bound in the image
    title(sprintf('Pick lower bound for the raster') )
    
    [xpos(4,1),ypos(4,1)]=ginput(1);
    hold on
    plot(xpos(1,1),ypos(1,1),'+r')
    
    R2=sqrt((xcross-xpos(4,1))*(xcross-xpos(4,1))+(ycross-ypos(4,1))*(ycross-ypos(4,1)))
    xcrossmatrix=ones(size(pos(:,selectedimage)))*xcross;
    ycrossmatrix=ones(size(pos(:,selectedimage)))*ycross;
    
    % Calculate Radius for all points
    Rall=sqrt((xcrossmatrix-pos(:,selectedimage)).*(xcrossmatrix-pos(:,selectedimage))+(ycrossmatrix-displ(:,selectedimage)).*(ycrossmatrix-displ(:,selectedimage)));
    
    % Calculate Angle for all points relativ to circle center
    newpos=pos(:,selectedimage)-xcross;newdispl=displ(:,selectedimage)-ycross;
    angleallpoints=acos((newpos.*xvector(1,1)+newdispl.*xvector(2,1))./(sqrt(newpos.*newpos+newdispl.*newdispl).*sqrt(xvector(1,1).*xvector(1,1)+xvector(2,1).*xvector(2,1))))*180/pi;
    negativangle=find(displ(:,selectedimage)<ycross);
    angleallpoints(negativangle)=angleallpoints(negativangle)*(-1)+360;
    selectpoints=find(Rall>min(R,R2) & Rall<max(R,R2) & angleallpoints>maxangle & angleallpoints<minangle);
    plot(pos(selectpoints,selectedimage),displ(selectpoints,selectedimage),'.r')
    drawnow
    
    
    % Do you want to keep the grid?
    confirmselectionraster = menu(sprintf('Do you want to use this raster?'),...
        'Yes','No, try again','Go back to raster-type selection');
    
    if confirmselectionraster==1
        rasterx=pos(selectpoints,:);
        rastery=displ(selectpoints,:);
        close all
    end
    
    if confirmselectionraster==2
        rasterx=0;
        rastery=0;
        close all
        hold off
        exist('rasterx')
        exist('rastery')
        exist('x')
        exist('y')
        exist('validx')
        exist('validy')
        [validx,validy,rasterx,rastery,selectedimage,x,y] = circ_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
    end
    
    if confirmselectionraster==3
        close all
        [rasterx,rastery,validx,validy,x,y]=gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);
    end
    
    
end    

% 
%-------------------------------
%

function [validx,validy,rasterx, rastery,selectedimage,x,y]=rect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);

[looppoints loopimages]=size(validx);
defaultimage=loopimages;

if x==1
    pos=validx;
end
if x==2
    pos=validy;
end
if x==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    pos=validx-validxfirst;
end
if x==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    pos=validy-validyfirst;
end
if y==1
    displ=validx;
end
if y==2
    displ=validy;
end
if y==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    displ=validx-validxfirst;
end
if y==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    displ=validy-validyfirst;
end

prompt = {'From which image do you want to select markers?'};
dlg_title = 'Marker selection';
num_lines= 1;
if selectedimage==0
    defaultimage=loopimages;
end
if selectedimage~0
    defaultimage=selectedimage;
end
def     = {num2str(defaultimage)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
selectedimage = str2num(cell2mat(answer(1,1)));
if selectedimage>loopimages
    selectedimage=loopimages;
end
if selectedimage<1
    selectedimage=1;
end

h=figure;
plot(pos(:,selectedimage),displ(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
title(sprintf('Define the region of interest.  Pick (single click) a point in the LOWER LEFT of your selection.\n  Do the same for a point in the UPPER RIGHT.'))
hold on
[xpos(1,1),ypos(1,1)]=ginput(1);
hold on
plot(xpos(1,1),ypos(1,1),'+b')
drawnow

[xpos(2,1),ypos(2,1)]=ginput(1);
hold on
plot(xpos(2,1),ypos(2,1),'+b')
drawnow

xmin = min(xpos);
xmax = max(xpos);
ymin = min(ypos);
ymax = max(ypos);

lowerline=[xmin ymin; xmax ymin];
upperline=[xmin ymax; xmax ymax];
leftline=[xmin ymin; xmin ymax];
rightline=[xmax ymin; xmax ymax];

plot(lowerline(:,1),lowerline(:,2),'-b')
plot(upperline(:,1),upperline(:,2),'-b')
plot(leftline(:,1),leftline(:,2),'-b')
plot(rightline(:,1),rightline(:,2),'-b')

selectpoints=find(pos(:,selectedimage)>min(xpos) & pos(:,selectedimage)<max(xpos) & displ(:,selectedimage)<max(ypos) & displ(:,selectedimage)>min(ypos))
clear pos
clear displ

if x==1
    pos=validx(selectpoints,:);
end
if x==2
    pos=validy(selectpoints,:);
end
if x==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    pos=validx(selectpoints,:)-validxfirst(selectpoints,:);
end
if x==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    pos=validy(selectpoints,:)-validyfirst(selectpoints,:);
end
if y==1
    displ=validx(selectpoints,:);
end
if y==2
    displ=validy(selectpoints,:);
end
if y==3
    validxfirst=zeros(size(validx));
    validxfirst=validx(:,1)*ones(1,loopimages);
    displ=validx(selectpoints,:)-validxfirst(selectpoints,:);
end
if y==4
    validyfirst=zeros(size(validy));
    validyfirst=validy(:,1)*ones(1,loopimages);
    displ=validy(selectpoints,:)-validyfirst(selectpoints,:);
end

plot(pos(:,selectedimage),displ(:,selectedimage),'o','MarkerEdgeColor','k','MarkerFaceColor','r')
title(sprintf('Red dots represent your new raster.'))
hold off

% Do you want to keep the grid?
confirmselection = menu(sprintf('Do you want to use this raster?'),...
    'Yes','No, try again','Go back to raster-type selection');

if confirmselection==1
    rasterx=pos;
    rastery=displ;
    close all
end

if confirmselection==2
    rasterx=0;
    rastery=0;
    close all
    hold off
    exist('rasterx')
    exist('rastery')
    exist('x')
    exist('y')
    exist('validx')
    exist('validy')
    [validx,validy,rasterx, rastery,selectedimage,x,y]=rect_grid_displ(validx,validy,x,y,selectedimage,rasterx,rastery);
end

if confirmselection==3
    close all
    [rasterx,rastery,validx,validy,x,y]=gridtypeselection(validx,validy,x,y,selectedimage,rasterx,rastery);
end
