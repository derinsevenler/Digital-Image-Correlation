% Script to save the displacement and strain matrices as images

load variables.mat

% Adjust the min and max of the strain map to match grid_generator
eles = dydx(find(isfinite(dydx)));
minim = min(min(eles));
maxim = max(max(eles));
range = maxim - minim;
scaled = (dydx-minim)*255/range;
scaled(~isfinite(scaled)) = 0;
scaled = uint8(scaled);
imwrite(scaled,'dydx.jpg','jpeg');

eles = dxdy(find(isfinite(dxdy)));
minim = min(min(eles));
maxim = max(max(eles));
range = maxim - minim;
scaled = (dxdy-minim)*255/range;
scaled(~isfinite(scaled)) = 0;
scaled = uint8(scaled);
imwrite(scaled,'dxdy.jpg','jpeg');


imwrite(im1cut,'im1cut.jpg','jpeg');


