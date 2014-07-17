clear all

[daq,PathStress] = uigetfile('*.txt','Open data aquisition file');
cd(PathStress);
import=importdata(daq,'\t');
clear daq
if exist('import.data')==0
    daq=import;
else
    daq=import.data;
end
prompt = {'Which column is the experimental time?','Which column is the experimental stress?','Which column is the experimental temperature?'};
dlg_title = 'Please specify the matrix properties';
num_lines= 1;
exptime  = num2str(2);
expstress  = num2str(7);
exptemp  = num2str(10);
def={exptime,expstress,exptemp};
answer = inputdlg(prompt,dlg_title,num_lines,def);
exptime = str2num(cell2mat(answer(1,1)));
expstress = str2num(cell2mat(answer(2,1)));
exptemp = str2num(cell2mat(answer(3,1)));

[Time_Image,PathImage] = uigetfile('*.txt','Open image-time file');
cd(PathImage);
Time_Image=importdata(Time_Image,'\t');
Time_Image(:,2)=Time_Image(:,2)-Time_Image(1,2)+1;

prompt = {'Time shift between data aquisition and image capture in [s]:'};
dlg_title = 'Please specify the onset of time';
num_lines= 1;
timeshift  = num2str(25);
def={timeshift};
answer = inputdlg(prompt,dlg_title,num_lines,def);
timeshift = str2num(cell2mat(answer(1,1)));

[Resultsfile,PathResults] = uiputfile('stress_image_x.txt','Where do you want to save the results?');

[loopmean with]=size(daq);

for i=5:(loopmean-5)
time_stress(i,1)=mean(daq(i-4:i,exptime));
time_stress(i,2)=mean(daq(i-4:i,expstress));
% time_stress(i,3)=mean(daq(i-4:i,exptemp));
end

[loopimage widthtime]=size(Time_Image);
for j=1:loopimage
    minpos=time_stress(:,1)-Time_Image(j,2)-timeshift;
    impos=find(abs(minpos)==min(abs(minpos)));
    image_stress(j,1)=j;
    image_stress(j,2)=mean(time_stress(impos,2));
    % image_stress(j,3)=mean(time_stress(impos,3));
end


cd(PathResults);
save  stress_image_x.txt image_stress -ASCII
plot(image_stress(:,1),image_stress(:,2))