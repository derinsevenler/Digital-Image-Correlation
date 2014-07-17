clear
close all
load variables.mat
load filenamelist.mat

to_the_left = true;

% Get location of strain map with relation to im1
offset = [min(min(YI)) min(min(XI))];
rscale = (max(max(YI))-min(min(YI)))/size(dydx,1);
cscale = (max(max(XI))-min(min(XI)))/size(dydx,2);

% Select the points of the gammaxy (dxdy) concentrations
factor = 1;
thresh = factor*std(dydx(isfinite(dydx)));
average = mean(dydx(isfinite(dydx)));
thresh2 = factor*std(dxdy(isfinite(dxdy)));
average2 = mean(dxdy(isfinite(dxdy)));
if to_the_left
    [slipr slipc] = find((dydx -average)>thresh);
    [slipr2 slipc2] = find((dxdy -average2)>thresh2);
else
    [slipr slipc] = find((dydx -average)<-thresh);
    [slipr2 slipc2] = find((dxdy -average2)<-thresh2);
end
slipr = slipr*rscale + offset(1);
slipc = slipc*cscale + offset(2);
slipr2 = slipr2*rscale + offset(1);
slipc2 = slipc2*cscale + offset(2);

imf=imread(filenamelist(end,:));
imshow(imf);
figure;
imshow(imf);
hold on;
plot(slipc, slipr,'*r');
title('Regions of Highest Shear strain \gamma_{xy}');

figure
imshow(imf);
hold on;
plot(slipc2, slipr2, '*r');
title('Regions of Highest Shear strain \gamma_{yx}');

