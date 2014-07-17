imagenum =0:54;

for i = 1:length(imagenum)
    imageNames{i} = strcat('inum', num2str(imagenum(i),'%04i'), '.tif');
end