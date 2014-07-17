%Written by Mark Buckley on 7/14/10
%DETERMINE AND SHOW STRAIN MAPS AFTER RUNNING DIC CODE

close all; clear all; 
scrsz = get(0,'ScreenSize');

%LOAD DATA
validx=dlmread('validx.dat');
validx(:,size(validx,2))=[];
validy=dlmread('validy.dat');
validy(:,size(validy,2))=[];

%SET UP displx
sizevalidx=size(validx);
validxfirst=zeros(size(validx));
validxfirst=mean(validx(:,1),2)*ones(1,sizevalidx(1,2));
displx=validx-validxfirst;

%SET UP disply
sizevalidy=size(validy);
validyfirst=zeros(size(validy));
validyfirst=mean(validy(:,1),2)*ones(1,sizevalidy(1,2));
disply=validy-validyfirst;

%MAKE MESH GRID

gridsizex=10*round(min(min(validx))/10):10:10*round(max(max(validx))/10);
gridsizey=10*round(min(min(validy))/10):10:10*round(max(max(validy))/10);
[XI,YI]=meshgrid(gridsizex,gridsizey);
gridstep=XI(1,2)-XI(1,1);


for last=2:1:length(validx(1,:))-1

    %CONFORM DISPLACEMENT DATA TO GRID
    Zx=griddata(validx(:,last),validy(:,last),displx(:,last),XI,YI,'cubic');
    Zy=griddata(validx(:,last),validy(:,last),disply(:,last),XI,YI,'cubic');

    %SET UP dydx, AKA gammaxy (5PLSQ)
    for i=1:1:length(Zy(:,1));
        dydx(i,:) = SlopeFinder(XI(i,:),Zy(i,:));
    end

    %SET UP dxdy AKA gammayx (5PLSQ)
    for j=1:1:length(Zx(1,:));
        dxdy(:,j) = SlopeFinder(YI(:,j),Zx(:,j));
    end

    cauchystrain=0.5*(dxdy+dydx);
    
    %SHOW RELEVANT IMAGE SECTION
    inum=81+last;
    iname=['inum',num2str(inum,'%04i'),'.tif'];
    I=imread(iname);
    Icut=I(min(min(YI)):max(max(YI)),min(min(XI)):max(max(XI)));
    figure('Color',[1 1 1],'Position', scrsz,'Windowstyle', 'docked');
    imagesc(Icut,[min(min(I)),max(max(I))]); colormap('gray'); axis off; 

    %GRAB FRAME FOR MOVIE
    flick(last-1)=getframe;
    close all;
    
    %PLOT HEAT MAPS
    if last==2;
        dydx1=dydx;
    end
    figure('Color',[1 1 1],'Position', scrsz,'Windowstyle', 'docked');
    imagesc(dydx, [min(min(dydx1)),max(max(dydx1))]); colormap('jet'); colorbar
    title('\gamma_{xy}','FontSize',30); axis off; set(colorbar,'fontsize',25);
    
    %GRAB FRAME FOR MOVIE    
    flick(last-1)=getframe;
    close all;
    
end

movie(flick,1);
movie2avi(flick,'montage.avi','compression','None','quality',100);
%movie2avi(flick,'strainmovie.avi','compression','Cinepak','quality',100); %'Indeo3','quality',50);
%movie2avi(flick,'strainmovie.avi','compression','None','quality',100);



