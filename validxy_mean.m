%% This function allows to take the mean of a couple of analyszed markers.
% Furthermore it allows you to take also the mean of the time. Every
% experimentalist knows the issue, vibrations. The only way of reducing the
% influence of vibrations which you have in your images is to take a mean
% value through time or images.

function [validx_mean,validy_mean,time_image_mean]=validxy_mean(validx,validy, time_image);


if exist('validx')==0
    [validxname,Pathvalidx] = uigetfile('*.dat','Open validx.dat');
    if validxname==0
        disp('You did not select a file!')
        return
    end
    cd(Pathvalidx);
    validx=importdata(validxname,'\t');
end
if exist('validy')==0
    [validyname,Pathvalidy] = uigetfile('*.dat','Open validy.dat');
    if validyname==0
        disp('You did not select a file!')
        return
    end
    cd(Pathvalidy);
    validy=importdata(validyname,'\t');
end
if exist('time_image')==0
    [timename,Pathtime] = uigetfile('*.txt','Open time_image.txt');
    if timename==0
        time=0;
    else
        cd(Pathtime);
        time_image=importdata(timename,'\t');
    end
end



[nummarker numimages]=size(validx);
numaverage=3;
for j=1:numaverage
    validxtemp(:,:,j)=validx(:,j:numaverage:(numimages-numaverage)+j);
    validytemp(:,:,j)=validy(:,j:numaverage:(numimages-numaverage)+j);
    if time~=0
        time_image_temp(:,:,j)=time_image(j:numaverage:(numimages-numaverage+1)+j,:);
    end
end
validx_mean=mean(validxtemp,3);
validy_mean=mean(validytemp,3);
save validx_mean.dat validx_mean
save validy_mean.dat validy_mean
if time~=0
    time_image_mean=mean(time_image_temp,3);
    save time_image.dat time_image
end