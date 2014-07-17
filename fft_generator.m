
close all;
load variables.mat;

field = dydx;
%imshow(im1);
%figure;
imagesc(field);
title('Select region: click lower-left corner and upper-right corner');
[x,y] = ginput(2);

c1 = min(x);
c2 = max(x);
r1 = min(y);
r2 = max(y);

fftim = field(r1:r2, c1:c2); 

% Adjust the range of fftim to go from 0 to 1
smallest = min(min(fftim));
largest = max(max(fftim));
fftim = (fftim+abs(smallest))./abs(largest-smallest); % fftim is of range 0-1.

% Subtract out mean to remove the DC imprint
fftim = fftim - mean(mean(fftim));

f = fftshift(fft2(fftim));


g = abs(f);
g = g/sum(g(:));
% g = uint8(g*255/max(max(g))); % Convert to uint8 grayscale.

[h,w] = size(g);
xs = linspace(-w/2,w/2,w);
ys = linspace(-h/2,h/2,h);

%f2 = figure;
%imshow(x,y,g, []); colormap('gray'); % Show the fft

f3 = figure;
[X,Y] = meshgrid(xs,ys);
contour(X,Y,flipud(g));
colorbar; % The matrix is flipped to maintain consistency with image coordinates
axis image
xlabel('cycles in image');
ylabel('cycles in image');

% Find the locations of the maximums in the strain map (there will be two, by
% symmetry)
peakval = max(max(g));
[peakr peakc] = find(g == peakval);

% Find the distance between the maximums and the angle of the line
% connecting them above the horizontal (radians)

deltax = peakc(2) - peakc(1);
deltay = peakr(2) - peakr(1);

dist = sqrt(deltax^2 + deltay^2);
angle = atan(deltay/deltax);

lambdax = w/deltax;
lambday = 9999;
if deltay
    lambday = h/deltay;
end

% TODO: Calculate 'sharpness' of peak
intensity = peakval./mean(mean(g));

lambdax
lambday
intensity