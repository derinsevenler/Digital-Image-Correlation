close all;


% Use strain_MRB calculated displacements
load variables.mat;
load grid_x.dat;
load grid_y.dat;

% Build a mesh of the true x and y displacements to match displx, disply

totalStrain = .09;
totaldisp = totalStrain*512;     % total displacement in pixels
m = totaldisp/512;              % displacement amount, per pixel from bottom
truedisplx = m*(mean(grid_y)*2- grid_y);
truedisply = zeros(size(grid_y));

% last (largest) displacement vector from displx, matches truedisplx coordinates
displacementx = displx(:,last);
displacementy = disply(:,last);
% flippedgrid_y = mean(grid_y)*2- grid_y;

positionx = grid_x + displacementx;
positiony = grid_y + displacementy;
truepositionx = grid_x + truedisplx;
truepositiony = grid_y + truedisply;

% Identify points that are outside a certain tolerance
tolerance = 5;  % say, at most 5 pixels off is acceptable

[badn] = find(sqrt((positionx - truepositionx).^2 +(positiony -truepositiony).^2) > tolerance);


figure; hold on;
plot(grid_x,grid_y, '+b');
plot(grid_x(badn),grid_y(badn),'*r');
legend('gridpoints within tolerance', 'Gridpoints outside tolerance');
hold off;

% Plot the expected displacement grid and the actual displacement grid

figure; hold on;
set(gca, 'YDir', 'reverse');
plot(positionx, positiony, 'ob');
plot(truepositionx, truepositiony, '*r');

hold off;

% Get all the positions and values in the strain map
[r,c] = find(~isnan(dxdy));
vals = dxdy(~isnan(dxdy));
[r2,c2] = find(~isnan(dydx));
vals2 = dydx(~isnan(dydx));

% Plot the strain measurements
figure;
plot(sqrt(r.^2 + c.^2),vals,'.');

title('Error Analysis - \gamma _y_x');
xlabel('gridpoint distance from upper left corner of image');
ylabel('\gamma _y_x');

figure;
plot(sqrt(r2.^2 + c2.^2),vals2,'.');

title('Error Analysis - \gamma _x_y');
xlabel('gridpoint distance from upper left corner of image');
ylabel('\gamma _x_y');

figure;
hist(vals,30);
title('\gamma _y_x spread');
ylabel('number');
xlabel('\gamma _y_x');

figure;
hist(vals2,30);
title('\gamma _x_y spread');
xlabel('\gamma _x_y');
ylabel('number');

% Build a 'good' gridpoint set, excluding  points outside tolerance

% % Import grid_x.dat, grid_y.dat
% load grid_x.dat;
% load grid_y.dat;
% % Save them as grid_x_old.dat, grid_y_old.dat
% save grid_x_old.dat grid_x -ascii -tabs
% save grid_y_old.dat grid_y -ascii -tabs
% 
% % build new gridpoints
% grid_x(badn)=[];
% grid_y(badn)=[];
% 
% % Save conditioned grid as grid_x.dat, grid_y.dat
% save grid_x.dat grid_x -ascii -tabs
% save grid_y.dat grid_y -ascii -tabs


% import the displacement  and position matrices:
% The coordinates of gridpoint n is stored in (grid_x(n),grid_y(n))
% The position of gridpoint n in image i is stored in
% (validx(n,i),validy(n,i)).

% % Compare original and strained image with original and strained true
% % points
% I = imread('inum0000.tif');
% figure; hold on;
% imagesc(I); colormap('gray'); axis off;
% plot(grid_x,grid_y,'+b');
% I = imread('inum0019.tif');
% figure; hold on;
% imagesc(I); colormap('gray'); axis off;
% plot(truepositionx,truepositiony,'+b');
