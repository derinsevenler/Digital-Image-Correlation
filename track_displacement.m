function G = track_displacement(area,thickness,microns_per_pixel,fps,hz,Newtons_per_micron,forcing_displacement,lastImage)
% area in mm^2, thickness and umresolution in um. imageNames is a cell array
% of strings, the sequential names of each image. hz is frequency in Hz.
close all;
thresh = 150; % threshold value for 'edge' to be identified

imagenum = 1:lastImage;
imageNames = {};
for i = 1:(lastImage)
    imageNames{i} = strcat('inum', num2str(imagenum(i),'%04i'), '.tif');
end

num = length(imageNames);
t = (1:num)'/fps; % timescale

% For each image, measure position of the edge, at the middle of the image
edges = zeros(num,1);
for i = 1:length(edges)
    I = imread(imageNames{i});
    section = I(250,:);
%     if i ==1
%         figure;
%         plot(section)
%     end
    edges(i) = find(section>thresh,1);
end

% Convert pixel measurements to microns
edges = (edges - min(edges)).*microns_per_pixel;

% calculate best-fit sinusoid with frequency hz

% General Linear Regression
X = [ones(num,1) cos((2*pi*hz)*(t)) sin((2*pi*hz)*(t))];
scoeffs = X\edges;

tfit = (1:.01:num)/fps;
tfit = tfit';
yfit = [ones(length(tfit),1) cos((2*pi*hz)*(tfit)) sin((2*pi*hz)*(tfit))]*scoeffs;

figure; hold on;
plot(t,edges, '*k');
plot(tfit, yfit,'-r');
xlabel('time'); ylabel('microns'); title('Sinusoidal Regression - Horizontal Displacement');

% Determine stiffness using section area, thickness, and force calibration
% Assuming homogeneous 
pkpk = range(yfit);
u = forcing_displacement-pkpk;  % in microns
gamma = u/thickness;        % microns/microns
F = Newtons_per_micron*pkpk;
tau = F/area*10^6;      % Newtons per m^2. 1 m^2 = 10^6 mm^2

G = tau/gamma;

end