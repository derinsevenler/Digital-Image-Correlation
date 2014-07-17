% Captures frames with your images and an overlay of the tracking marker
% written by Chris

function [validx,validy]=Markerplotting(validx,validy);

if exist('validx')==0
    [validx,Pathvalidx] = uigetfile('*.mat','Open validx.mat');
    cd(Pathvalidx);
    validx=importdata(validx,'\t');
end
if exist('validy')==0
    [validy,Pathvalidy] = uigetfile('*.mat','Open validy.mat');
    cd(Pathvalidy);
    validy=importdata(validy,'\t');
end
if exist('filenamelist')==0
load('filenamelist')
end

sizevalidx=size(validx);
sizevalidy=size(validy);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);

h=figure;
videoselection = menu(sprintf('Do you want to create a video?'),'Yes','No')
if videoselection==1
    mkdir('Video_Markers')
    Vid='Vid';
end
for i=1:(loopimages-1)
    imshow(filenamelist(i,:))
    hold on
    title(['Displacement in y-direction versus x-y-position',sprintf(' (Current image #: %1g)',i)]);
    plot(validx(:,i),validy(:,i),'.g','MarkerSize',2)
    hold off
    drawnow
    if videoselection==1
        u=i+10000;
        ustr=num2str(u);
        videoname=[Vid ustr '.jpg']
        cd('Video_Markers');
        saveas(h,videoname,'jpg')
        cd('..')
    end
    
end