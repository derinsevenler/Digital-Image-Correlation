%Written by Mark Buckley on 7/14/10
%DETERMINE AND SHOW STRAIN MAPS AFTER RUNNING DIC CODE
function strain_DDSa

close all; clear all; 

%LOAD DATA

if exist('validx_corr.dat','file')
    validx=dlmread('validx_corr.dat');
else
    validx=dlmread('validx.dat');
end
validx(:,size(validx,2))=[];
if exist('validy_corr.dat','file')
    validy=dlmread('validy_corr.dat');
else
    validy=dlmread('validy.dat');
end
validy(:,size(validy,2))=[];

load filenamelist.mat

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
last=length(validx(1,:));
%gridsizex=10*round(min(min(validx))/10):10:10*round(max(max(validx))/10);
%gridsizey=10*round(min(min(validy))/10):10:10*round(max(max(validy))/10);
gridsizex=10*round(min(validx(:,last))/10):5:10*round(max(validx(:,last))/10);
gridsizey=10*round(min(validy(:,last))/10):5:10*round(max(validy(:,last))/10);
[XI,YI]=meshgrid(gridsizex,gridsizey);
gridstep=XI(1,2)-XI(1,1);

%CONFORM DISPLACEMENT DATA TO GRID
Zx=griddata(validx(:,last),validy(:,last),displx(:,last),XI,YI,'cubic');
Zy=griddata(validx(:,last),validy(:,last),disply(:,last),XI,YI,'cubic');


% %PLOT DISPLACEMENT MESH
% figure; mesh(XI,YI,Zx); 
% figure; mesh(XI,YI,Zy); 

% %SET UP dydx (FORWARD DIFFERENCE)
% for i=1:1:length(Zy(:,1));
%     dydx(i,length(Zy(1,:)))=NaN;
%     for j=1:1:length(Zy(1,:))-1;
%         dydx(i,j)=Zy(i,j+1)-Zy(i,j);
%     end;        
% end

% %SET UP dxdy (FORWARD DIFFERENCE)
% for j=1:1:length(Zx(1,:));
%     dxdy(length(Zx(:,1)),j)=NaN;
%     for i=1:1:length(Zx(:,1))-1;
%         dxdy(i,j)=Zx(i+1,j)-Zx(i,j);
%     end;
% end

%SET UP dydx, AKA gammaxy (5PLSQ) 
for i=1:1:length(Zy(:,1));
    dydx(i,:) = SlopeFindera(XI(i,:),Zy(i,:)); 
end

%SET UP dxdy AKA gammayx (5PLSQ)
for j=1:1:length(Zx(1,:));
    dxdy(:,j) = SlopeFindera(YI(:,j),Zx(:,j));    
end
%dydx = fix_outliers(dydx);
%dxdy = fix_outliers(dxdy);

cauchystrain=0.5*(dxdy+dydx);

% %PLOT STRAIN MESH
% 
% figure; mesh(XI,YI,dxdy); 
% figure; mesh(XI,YI,dydx); 

%LOAD AND PLOT THE RELEVANT SECTION OF THE IMAGE
%XI(:,1) = ones;
im1=imread(filenamelist(1,:));
im1cut=im1(min(min(YI)):max(max(YI)),min(min(XI)):max(max(XI)));
imf=imread(filenamelist(end,:));
imfcut=imf(min(min(YI)):max(max(YI)),min(min(XI)):max(max(XI)));
figure('Color',[1 1 1]);imagesc(im1cut); colormap('gray'); axis off;
%    plot(validxfirst,validyfirst,'g+')  %OVERLAY GRID 
%figure('Color',[1 1 1]);imagesc(imfcut,[min(min(imfcut)),max(max(imfcut))]); colormap('gray'); axis off;
%    plot(validxfirst,validyfirst,'g+')  %OVERLAY GRID 

%PLOT HEAT MAPS

xmicrons = (max(max(XI)) - min(min(XI))) * 1.3;
ymicrons = (max(max(YI)) - min(min(YI))) * 1.3;
x = linspace(0,xmicrons,size(Zx,2));
y = linspace(0,ymicrons,size(Zy,2));

figure('Color',[1 1 1]); imagesc(x,y,Zx,[min(min(Zx)),max(max(Zx))]); colormap('jet'); colorbar
title('u_{x}','FontSize',30); axis image; set(colorbar,'fontsize',18);
figure('Color',[1 1 1]); imagesc(x,y,Zy,[min(min(Zy)),max(max(Zy))]); colormap('jet'); colorbar
title('u_{y}','FontSize',30);axis image; set(colorbar,'fontsize',18);
figure('Color',[1 1 1]); imagesc(x,y,dxdy, [min(min(dxdy)),max(max(dxdy))]); colormap('jet'); colorbar
title('\gamma_{yx}','FontSize',30); axis image; set(colorbar,'fontsize',18);
figure('Color',[1 1 1]); imagesc(x,y,dydx, [min(min(dydx)),max(max(dydx))]); colormap('jet'); colorbar
title('\gamma_{xy}','FontSize',30); axis image; set(colorbar,'fontsize',18);
xlabel('\mum')
ylabel('\mum')

save variables.mat
end
