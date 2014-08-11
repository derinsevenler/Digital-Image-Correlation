%Written by Mark Buckley on 7/14/10
%DETERMINE AND SHOW STRAIN MAPS AFTER RUNNING DIC CODE
function strain_DDS

close all; clear all; 

%LOAD DATA

load grid_x.dat;
load grid_y.dat;

if exist('validx_corr.dat','file')
    validx=dlmread('validx_corr.dat');
else
    validx=dlmread('validx.dat');
end

if exist('validy_corr.dat','file')
    validy=dlmread('validy_corr.dat');
else
    validy=dlmread('validy.dat');
end


load filenamelist.mat;

%SET UP displx, disply
num_cp=size(validx,1);
displx=validx(:,end)-validx(:,1);
disply=validy(:,end)-validy(:,1);

% Make the Mesh
gridsizex=min(validx(:,end)):max(validx(:,end));
gridsizey=min(validy(:,end)):max(validy(:,end));
[XI,YI]=meshgrid(gridsizex,gridsizey);
gridstep=XI(1,2)-XI(1,1);

%CONFORM DISPLACEMENT DATA TO GRID
Zx=griddata(validx(:,end),validy(:,end),displx(:,end),XI,YI,'cubic');
Zy=griddata(validx(:,end),validy(:,end),disply(:,end),XI,YI,'cubic');

Zx = fix_outliers(Zx);
Zy = fix_outliers(Zy);

dxdx = zeros(size(Zx));
dxdy = zeros(size(Zx));
dydx = zeros(size(Zx));
dydy = zeros(size(Zx));

% Calculate normal and shear strains
for i=1:length(Zy(:,1))
    dxdx(i,:) = SlopeFindera(XI(i,:),Zx(i,:));    
    dydx(i,:) = SlopeFindera(XI(i,:),Zy(i,:));
end

for j=1:length(Zx(1,:));
    dydy(:,j) = SlopeFindera(YI(:,j),Zy(:,j));
    dxdy(:,j) = SlopeFindera(YI(:,j),Zx(:,j)); 
end

save variables.mat

% Plot Figures
figure('Color',[1 1 1]); imagesc(Zx); colormap('jet'); colorbar
title('u_{x}','FontSize',30); axis off;
figure('Color',[1 1 1]); imagesc(Zy); colormap('jet'); colorbar
title('u_{y}','FontSize',30);axis off;

figure('Color',[1 1 1]); imagesc(dxdx); colormap('jet'); colorbar
title('\epsilon _x','FontSize',16); axis off;
figure('Color',[1 1 1]); imagesc(dydy); colormap('jet'); colorbar
title('\epsilon _y','FontSize',16); axis off;
figure('Color',[1 1 1]); imagesc(dydx); colormap('jet'); colorbar
title('\gamma_x_y','FontSize',16); axis off;
figure('Color',[1 1 1]); imagesc(dxdy); colormap('jet'); colorbar
title('\gamma_y_x','FontSize',16); axis off;

save variables.mat
end
