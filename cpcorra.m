function xyainput = cpcorra(varargin)
%CPCORR Tune control point locations using cross-correlation. 
%   INPUT_POINTS = CPCORR(INPUT_POINTS_IN,BASE_POINTS_IN,INPUT,BASE) uses
%   normalized cross-correlation to adjust each pair of control points
%   specified in INPUT_POINTS_IN and BASE_POINTS_IN.
%
%   INPUT_POINTS_IN must be an M-by-3 double matrix containing the
%   coordinates of control points in the input image.  BASE_POINTS_IN is
%   an M-by-2 double matrix containing the coordinates of control points
%   in the base image. Third col of input_points_in is angle in radians
%
%   CPCORR returns the adjusted control points in INPUT_POINTS, a double
%   matrix the same size as INPUT_POINTS_IN.  If CPCORR cannot correlate a
%   pairs of control points, INPUT_POINTS will contain the same coordinates
%   as INPUT_POINTS_IN for that pair.

%   INPUT:   3-D, real, full matrix
%            logical, uint8, uint16, or double
%            must be finite (no NaNs, no Infs inside regions being
%            correlated). a_input_in is the rotation angles array, degrees.
%
%   BASE:    2-D, real, full matrix
%            logical, uint8, uint16, or double
%            must be finite (no NaNs, no Infs inside regions being
%            correlated)
%
%   For more information, see >>help cpcorr

[xyinput_in,a_input_in,xybase_in,input,base,last_input,doRot] = ParseInputs(varargin{:});

CORRSIZE = 10;
if doRot
    ANG_CHECK = [0 -2 2 -4 4 -6 6];
else
    ANG_CHECK = 0;
end

% get the intial large rectangles, centered at xyinput_in
% which are twice as large as what we'll use to do the actual image
% correlation - we'll crop to this size, then rotate, then crop again.
if doRot
    rects_input = calc_rects(xyinput_in,2.0*CORRSIZE,input);
else
    rects_input = calc_rects(xyinput_in,1.5*CORRSIZE,input);
end
% Get rectangles for the base image - these don't get rotated.
rects_base = calc_rects(xybase_in,CORRSIZE,base);
% Get rectangles for last input image
rects_last_input = calc_rects(xyinput_in, CORRSIZE,last_input);

ncp = size(xyinput_in,1);

xyinput = xyinput_in; % initialize adjusted control points matrix
a_input = a_input_in;
xyainput = [xyinput a_input];

for icp = 1:ncp
    
    % Get the image section of base which we want to find  'sub_image' in.
    sub_base_orig = imcrop(base,rects_base(icp,:));
    sub_last_input = imcrop(last_input,rects_last_input(icp,:));
    

    sub_base = (4*sub_base_orig + 0*sub_last_input)/4;
    
    topamp = 0;
    angle_best = a_input_in(icp);
    xpeak_best = 0;
    ypeak_best = 0;
    for rotangle = ANG_CHECK
        if doRot
            sub_input = imcrop(input,rects_input(icp,:));
            sub_input = imrotate(sub_input, rotangle);
            center = size(sub_input,1)/2;
            sub_input = imcrop(sub_input,[center-CORRSIZE*1.5, center-CORRSIZE*1.5, 3*CORRSIZE, 3*CORRSIZE]);

        else
            sub_input = imcrop(input,rects_input(icp,:));
        end
%         inputsize = size(sub_input);
%         basesize = size(sub_base);

        % The magical step - the normalized 2D cross-correlation
        norm_cross_corr = normxcorr2(sub_base,sub_input);
        s = size(norm_cross_corr);
        y2 = floor(CORRSIZE*2.5);
        x2 = y2; % For square elements, which we are using. 
        
        % find the peak of this surface to subpixel resolution
        [xpeak, ypeak, amplitude] = findpeak(norm_cross_corr,true);
        if amplitude>topamp
            topamp = amplitude;
            angle_best = rotangle;
            xpeak_best = xpeak;
            ypeak_best = ypeak;
        end
    end
    
    if doRot
        theta = (a_input_in(icp) + angle_best)*pi/180;
        n = size(sub_input,1);
        if theta > 0
            xpeak_best = xpeak_best*cos(-theta) + ypeak_best*sin(-theta);
            ypeak_best = ypeak_best*cos(-theta) + (n-xpeak_best)*sin(-theta);
        elseif theta < 0
            xpeak_best = xpeak_best*cos(-theta) + (n-ypeak_best)*sin(-theta);
            ypeak_best = ypeak_best*cos(-theta) - xpeak_best*sin(-theta);
        end
    end

    % eliminate any poor correlations
%     THRESHOLD = 0.2;
%     if (topamp < THRESHOLD) 
%         % low correlation, unable to adjust
%         continue
%     end
    
    % offset found by cross correlation
    inputsize = size(sub_input);
    corr_offset = [ (xpeak_best-x2) (ypeak_best-y2) ];
    
    % eliminate any big changes in control points
    ind = find(abs(corr_offset) > (CORRSIZE-1));
    if ~isempty(ind)
        % peak of norxcorr2 not well constrained, unable to adjust
%         disp('large displacement');
%         icp
        continue
    end
    input_fractional_offset = xyinput(icp,:) - round(xyinput(icp,:)*100)/100;
    base_fractional_offset = xybase_in(icp,:) - round(xybase_in(icp,:)*100)/100;    
    
    % adjust control point
    xyinput(icp,:) = xyinput(icp,:) + corr_offset;
    a_input(icp) = a_input(icp) + angle_best;
    
    xyainput = [xyinput a_input];
end

%-------------------------------
%
function rect = calc_rects(xy,halfwidth,img)

% Calculate rectangles so imcrop will return image with xy coordinate inside center pixel

default_width = 2*halfwidth;
default_height = default_width;

% xy specifies center of rectangle, need upper left
upperleft = round(xy) - halfwidth;

% need to modify for pixels near edge of images
upper = upperleft(:,2);
left = upperleft(:,1);
lower = upper + default_height;
right = left + default_width;
width = default_width * ones(size(upper));
height = default_height * ones(size(upper));

% check edges for coordinates outside image
[upper,height] = adjust_lo_edge(upper,1,height);
[dum,height] = adjust_hi_edge(lower,size(img,1),height);
[left,width] = adjust_lo_edge(left,1,width);
[dum,width] = adjust_hi_edge(right,size(img,2),width);

% set width and height to zero when less than default size
iw = find(width<default_width);
ih = find(height<default_height);
idx = unique([iw; ih]);
width(idx) = 0;
height(idx) = 0;

rect = [left upper width height];

%-------------------------------
%
function [coordinates, breadth] = adjust_lo_edge(coordinates,edge,breadth)

indx = find( coordinates<edge );
if ~isempty(indx)
    breadth(indx) = breadth(indx) - abs(coordinates(indx)-edge);
    coordinates(indx) = edge;
end

%-------------------------------
%
function [coordinates, breadth] = adjust_hi_edge(coordinates,edge,breadth)

indx = find( coordinates>edge );
if ~isempty(indx)
    breadth(indx) = breadth(indx) - abs(coordinates(indx)-edge);
    coordinates(indx) = edge;
end

%-------------------------------
%
function [xyinput_in,a_input_in,xybase_in,input,base,last_input,doRot] = ParseInputs(varargin)

iptchecknargin(6,6,nargin,mfilename);

xyainput_in = varargin{1};
xyinput_in = xyainput_in(:,1:2);
a_input_in = xyainput_in(:,3);
xybase_in = varargin{2};
if size(xyinput_in,2) ~= 2 || size(xybase_in,2) ~= 2
    msg = sprintf('input must be M-by-3, base must be M-by-2');
    eid = sprintf('Images:%s:cpMatrixMustBeMby3',mfilename);
    error(eid,msg);
end

if size(xyinput_in,1) ~= size(xybase_in,1)
    msg = sprintf('In function %s, INPUT and BASE images need same number of control points.',mfilename);
    eid = sprintf('Images:%s:needSameNumOfControlPoints',mfilename);    
    error(eid,msg);
end

input = varargin{3};
base = varargin{4};
last_input = varargin{5};
doRot = varargin{6};

if ndims(input) ~= 2 || ndims(base) ~= 2
    msg = sprintf('In function %s, Images must be intensity images.',mfilename);
    eid = sprintf('Images:%s:intensityImagesReq',mfilename);        
    error(eid,msg);
end

input = double(input);
base = double(base);
last_input = double(last_input);

if any(xyinput_in(:,1)<0.5) || any(xyinput_in(:,2)<0.5) || any(xyinput_in(:,1)>size(input,2)+0.5) || ...
   any(xyinput_in(:,2)>size(input,1)+0.5) || ...
   any(xybase_in(:,1)<0.5) || any(xybase_in(:,2)<0.5) || any(xybase_in(:,1)>size(base,2)+0.5) || ...
   any(xybase_in(:,2)>size(base,1)+0.5)
    msg = sprintf('In function %s, Control Points must be in pixel coordinates.',mfilename);
    eid = sprintf('Images:%s:cpPointsMustBeInPixCoord',mfilename);
    error(eid,msg);
end
