% Initialize variables

n = 20; % Number of images to use
etot = .10; % Amount of shear strain to use

I = imread('image050.tif');
[h,w] = size(I);

step = 0;
for estep = linspace(0,etot,n)
    
    % Create sheared image
    T = maketform('affine',[1 0 0; -estep 1 0; 0 0 1]);
    R = makeresampler({'cubic','nearest'},'fill');
    B = imtransform(I,T,R);

    % Crop left side of shear image away to normalize resolution
     B = B(:, 1:w);
    % Or right side
%     if estep >0
%         B = B(:,(size(B,2)-w:size(B,2)));
%     end
    
    % Save sheared image as a numbered filename
    num = sprintf('%04.0f',step);
    filename = strcat('inum', num,'.tif');
    imwrite(B,filename,'TIFF');
    step = step+1;
end