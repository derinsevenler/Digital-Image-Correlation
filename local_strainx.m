%% This function calculates the local resolved strain in the horizontal or
% x-direction. It is recommended to run displacement.m first to clean up
% your validx.dat and validy.dat. local_displacement.m allows you to
% define e.g. 20 regions along x where the function calculates the gradient
% (strainx) and plots it for each image. The gradientx e.g. contains the
% strainx in the columns for each image in rows.

function [marker_displacementx marker_displacementy gradientx gradienty]=...
    local_strainx(validx, validy)

if exist('validx')==0
    [validx,Pathvalidx] = uigetfile('*.dat','Open validx.dat');
    cd(Pathvalidx);
    validx=importdata(validx,'\t');
end
if exist('validy')==0
    [validy,Pathvalidy] = uigetfile('*.dat','Open validy.dat');
    cd(Pathvalidy);
    validy=importdata(validy,'\t');
end

sizevalidx=size(validx);
sizevalidy=size(validy);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);

displacementx=(mean(validx(:,1))-mean(validx));
displacementy=(mean(validy(:,1))-mean(validy));


prompt = {'From which image do you want to select markers?'};
dlg_title = 'Marker selection';
num_lines= 1;
defaultimage=1;
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

plot(validx(:,selectedimage),validy(:,selectedimage),'+b');
title(sprintf('Define the region of interest.  Pick (single click) a point in the LOWER LEFT region of the gage section.\n  Do the same for a point in the UPPER RIGHT portion of the gage section.'))
hold on
[x(1,1),y(1,1)]=ginput(1);
plot(x(1,1),y(1,1),'+b')
plot([min(validx(:,selectedimage)); max(validx(:,selectedimage))], [y(1,1); y(1,1)],'-r')
plot([x(1,1),x(1,1)], [min(validy(:,selectedimage)),max(validy(:,selectedimage))],'-r')

[x(2,1),y(2,1)]=ginput(1);
hold on
plot(x(2,1),y(2,1),'+b')
plot([min(validx(:,selectedimage)); max(validx(:,selectedimage))], [y(2,1); y(2,1)],'-r')
plot([x(2,1),x(2,1)], [min(validy(:,selectedimage)),max(validy(:,selectedimage))],'-r')

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

lowerline=[xmin ymin; xmax ymin];
upperline=[xmin ymax; xmax ymax];
leftline=[xmin ymin; xmin ymax];
rightline=[xmax ymin; xmax ymax];

plot(lowerline(:,1),lowerline(:,2),'-g')
plot(upperline(:,1),upperline(:,2),'-g')
plot(leftline(:,1),leftline(:,2),'-g')
plot(rightline(:,1),rightline(:,2),'-g')

selectpoints=find(validx(:,selectedimage)>min(x) & validx(:,selectedimage)<max(x) & validy(:,selectedimage)<max(y) & validy(:,selectedimage)>min(y));

validx_new=validx(selectpoints,:);
validy_new=validy(selectpoints,:);

hold on

plot(validx_new(:,selectedimage),validy_new(:,selectedimage),'+g')
title(sprintf('Red dots represent your new raster.'))
hold off
drawnow


% ---------------------------------
gridselection1 = menu(sprintf('Do you want to use the green highlighted markers?'),...
    'Yes','No');

if gridselection1==2
    return
end

%% here you can input into how many areas you would like the strain to be
% split up

prompt = {'Into how many areas do you want to split the markers?'};
dlg_title = 'Divide markers';
num_lines= 1;
defaultsplit=10;
def     = {num2str(defaultsplit)};
answer = inputdlg(prompt,dlg_title,num_lines,def);
selectesplit = str2num(cell2mat(answer(1,1)));


if selectesplit<1
    selectesplit=20;
end
clear marker_displacement

hold on
for i=1:selectesplit
    clear validx_plot
    clear validy_plot
    clear selectmarkers
    posx(i)=xmin+((xmax-xmin)/selectesplit)*i;
    plot([posx(i);posx(i)], [ymin; ymax],'-m')
    selectmarkers=find(validx_new(:,selectedimage)>(xmin+((xmax-xmin)/selectesplit)*(i-1)) & validx_new(:,selectedimage)<(xmin+((xmax-xmin)/selectesplit)*(i)));
    validx_plot=validx_new(selectmarkers,:);
    validy_plot=validy_new(selectmarkers,:);
    plot(validx_plot(:,selectedimage),validy_plot(:,selectedimage),'+r')
    drawnow
    clear validy_local
    validy_local=validy_plot;
    marker_displacementy(:,i)=mean(validy_local)';
    clear validx_local
    validx_local=validx_plot;
    marker_displacementx(:,i)=mean(validx_local)';
    plot(validx_local(:,selectedimage),validy_local(:,selectedimage),'xg')
    drawnow
    drawnow
end

hold off

%% save everything

save marker_displacementx.txt marker_displacementx -ASCII
save marker_displacementy.txt marker_displacementy -ASCII

sizemvalidx=size(marker_displacementx);
clear mdisplx;
clear mdisply;
mvalidxfirst=ones(sizemvalidx(1,1),1)*marker_displacementx(1,:);
mdisplx=marker_displacementx-mvalidxfirst;
mvalidyfirst=ones(sizemvalidx(1,1),1)*marker_displacementy(1,:);
mdisply=marker_displacementy-mvalidyfirst;
clear mvalidxfirst
clear mvalidyfirst

grad_mdisplx=gradient(mdisplx);
grad_mdisply=gradient(mdisply);

grad_marker_displacementx=gradient(marker_displacementx);
grad_marker_displacementy=gradient(marker_displacementy);

gradientx=grad_mdisplx./grad_marker_displacementx;
gradienty=grad_mdisply./grad_marker_displacementy;

save gradientx.txt gradientx -ASCII
save gradienty.txt gradienty -ASCII

%% Here you can let the function plot the displacement

checkdisplacement = menu(sprintf('Do you want to check the displacement x?'),...
    'Yes','No');

if checkdisplacement==1


    minmeanx=min(min(marker_displacementx));
    maxmeanx=max(max(marker_displacementx));
    minmdisplx=min(min(mdisplx));
    maxmdisplx=max(max(mdisplx));

    figure
    for i=1:loopimages
        plot(marker_displacementx(i,:),mdisplx(i,:),'.r')
        title(['Displacement',sprintf(' (Current image #: %1g)',i)])
        axis([minmeanx maxmeanx minmdisplx maxmdisplx])
        drawnow
    end
end


%% Here, you can let the function plot the strainx
checkdisplacement = menu(sprintf('Do you want to check the local strain x?'),...
    'Yes','No');
videoselection = menu(sprintf('Do you want to create a video?'),'Yes','No');
if videoselection==1
    close all
    mkdir('video_local_strain')
    cd('video_local_strain');
    Vid='Vid';
end
if checkdisplacement==1
    minmeanx=min(min(marker_displacementx));
    maxmeanx=max(max(marker_displacementx));
    min_marker_gradx=min(min(gradientx));
    max_marker_gradx=max(max(gradientx));
    H=figure;
    for i=1:loopimages
        plot(marker_displacementx(i,:),gradientx(i,:),'.-b')
        title(['Gradient of Displacement',sprintf(' (Current image #: %1g)',i)])
        axis([minmeanx maxmeanx min_marker_gradx max_marker_gradx])
        xlabel('Position in Pixel [ ]')
        ylabel('Local Displacement Gradient [ ]')
        drawnow
        
        if videoselection==1
            u=i+10000;
            ustr=num2str(u);
            videoname=[Vid ustr '.jpg'];
            saveas(h,videoname,'jpg');
        end
        %                 pause(0.1)
    end
    if videoselection==1
        cd('..')
    end
end

