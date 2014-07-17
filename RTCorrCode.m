function [validx, validy, displx, disply]=RTCorrCode(grid_x,grid_y,straindir,Firstimagename)

% Real time Correlation Code
%
% Written by Chris

RTselection = menu(sprintf('End processing by end.txt or by last image?'),...
    'Stop with end.txt','Stop with image check','Exit');

if RTselection==1
end

if RTselection==2
end

if RTselection==3
    return
end


% Filename

if exist('Firstimagename')==0
    [Firstimagename ImageFolder]=uigetfile('*.tif','Open First Image');
    if Firstimagename~~[]
        cd(ImageFolder);
    end
end

if Firstimagename~~[]
% Get the number of image name
letters=isletter(Firstimagename);
Pointposition=findstr(Firstimagename,'.');
Firstimagenamesize=size(Firstimagename);
counter=Pointposition-1;
counterpos=1;
letterstest=0;
while letterstest==0
    letterstest=letters(counter);
    if letterstest==1
        break
    end
    Numberpos(counterpos)=counter;
    counter=counter-1;
    counterpos=counterpos+1;
    if counter==0
        break
    end
end

Filename_first = Firstimagename(1:min(Numberpos)-1);
Firstfilenumber=Firstimagename(min(Numberpos):max(Numberpos));
Lastname_first = Firstimagename(max(Numberpos)+1:Firstimagenamesize(1,2));
Firstfilenumbersize=size(Firstfilenumber);
onemore=10^(Firstfilenumbersize(1,2));
filenamelist(1,:)=Firstimagename;
h=figure;
if exist('grid_x')==0
    fpstest=1;
    Filelist=[Firstimagename;Firstimagename];
    while fpstest==1
        [grid_x,grid_y]=grid_generator(Firstimagename,ImageFolder);
        [processingtime]=fpstestfunc(grid_x,grid_y,Filelist);
        fpstest = menu(sprintf(['Processing the selected grid will allow ' , num2str(1/processingtime),' frames per second' ]),'Try again','Use the grid');
        if fpstest==1
            clear grid_x; clear grid_y;
        end
    end
end

Firstfilenumber=str2num(Firstfilenumber);
u=1+onemore+Firstfilenumber;
ustr=num2str(u);
filenamelist(2,:)=[Filename_first ustr(2:Firstfilenumbersize(1,2)+1) Lastname_first];
numberofimages=2;

counter=1;

input_points_x=grid_x;
input_points_y=grid_y;
base_points_x=grid_x;
base_points_y=grid_y;
base = uint8(mean(double(imread(filenamelist(1,:))),3));            % read in the base image ( which is always  image number one. You might want to change that to improve correlation results in case the light conditions are changing during the experiment
numberofmarkers=max(size(grid_x))*min(size(grid_x));
validx(:,1)=reshape(grid_x,[],1);
displx=zeros(numberofmarkers,1);
validy(:,1)=reshape(grid_y,[],1);
disply=zeros(numberofmarkers,1);
tic

while exist('end.txt','file') ==0;
    pause(0.01);

    if exist(filenamelist((counter+1),:),'file') ==2;
        warning(['# Processed Images: ', num2str(numberofimages-1),'; # markers:',num2str(numberofmarkers), '; Processing Image: ',filenamelist(counter+1,:)])    % plot a title onto the image

        input = mean(double(imread(filenamelist((counter+1),:))),3);       % read in the image which has to be correlated

        input_points_for(:,1)=reshape(input_points_x,[],1);         % we reshape the input points to one row of values since this is the shape cpcorr will accept
        input_points_for(:,2)=reshape(input_points_y,[],1);
        base_points_for(:,1)=reshape(base_points_x,[],1);
        base_points_for(:,2)=reshape(base_points_y,[],1);
        input_correl(:,:)=cpcorr(round(input_points_for), round(base_points_for), input, base);           % here we go and give all the markers and images to process to cpcorr.m which ic a function provided by the matlab image processing toolbox
        input_correl_x=input_correl(:,1);                                       % the results we get from cpcorr for the x-direction
        input_correl_y=input_correl(:,2);                                       % the results we get from cpcorr for the y-direction

        validx(:,counter+1)=input_correl_x;                                                     % lets save the data
        savelinex=input_correl_x';
        dlmwrite('resultsimcorrx.txt', savelinex , 'delimiter', '\t', '-append');       % Here we save the result from each image; if you are desperately want to run this function with e.g. matlab 6.5 then you should comment this line out. If you do that the data will be saved at the end of the correlation step - good luck ;-)

        validy(:,counter+1)=input_correl_y;
        saveliney=input_correl_y';
        dlmwrite('resultsimcorry.txt', saveliney , 'delimiter', '\t', '-append');

        base_points_x=grid_x;
        base_points_y=grid_y;
        input_points_x=input_correl_x;
        input_points_y=input_correl_y;

        subplot(2,2,1)
        imshow(filenamelist(counter+1,:))                     % update image
        hold on
        plot(grid_x,grid_y,'g+')                                % plot start position of raster
        plot(input_correl_x,input_correl_y,'r+')        % plot actual postition of raster
        hold off
        drawnow

        displx(:,counter+1)=validx(:,counter+1)-validx(:,1);
        disply(:,counter+1)=validy(:,counter+1)-validy(:,1);

        subplot(2,2,2)
        xdata=validx(:,counter+1);
        ydata=displx(:,counter+1);
        if counter==1
            x(1)=0
            x(2)=0
        end
        [x,resnormx,residual,exitflagx,output]  = lsqcurvefit(@linearfit, [x(1) x(2)], xdata, ydata);
        plot(xdata,ydata,'.');
        hold on;
        ydatafit=x(1)*xdata+x(2);
        plot(xdata,ydatafit,'r');
        hold off
        xlabel('x-pos [pixel]')
        ylabel('x-displ [pixel]')
        title('x displ. versus x pos. in [pixel]')

        slopex(counter,:)=[i x(1) x(2)];

        subplot(2,2,4)
        xdata=validy(:,counter+1);
        ydata=disply(:,counter+1);
        if counter==1
            y(1)=0
            y(2)=0
        end
        [y,resnormx,residual,exitflagx,output]  = lsqcurvefit(@linearfit, [y(1) y(2)], xdata, ydata);
        plot(xdata,ydata,'.g');
        hold on;
        ydatafit=y(1)*xdata+y(2);
        plot(xdata,ydatafit,'r');
        hold off
        xlabel('y-pos [pixel]')
        ylabel('y-displ [pixel]')
        title('y displ. versus y pos. in [pixel]')

        slopey(counter,:)=[i y(1) y(2)];

        subplot(2,2,3)
        plot(slopex(:,2),'-b')
        hold on
        plot(slopey(:,2),'-g')
        hold off
        xlabel('Image # [ ]')
        ylabel('x- and y-strain [ ]')
        title('Strain in x and y direction versus Image #')

        counter=counter+1;

        u=1+u;
        ustr=num2str(u);
        filenamelist(counter+1,:)=[Filename_first ustr(2:Firstfilenumbersize(1,2)+1) Lastname_first];
        [numberofmarkers numberofimages]=size(validx);
        
        if RTselection==2
            if exist(filenamelist((counter+1),:),'file') ==0;
                save validx.dat validx -ascii -tabs
                save validy.dat validy -ascii -tabs
                warning('Last image detected, RTCorrCode stopped')
                return
            end
        end
        
        
        subplot(2,2,1),title(['# Processed Images: ', num2str(numberofimages-1),'; fps: ', num2str((numberofimages-1)/toc),'; # markers:',num2str(numberofmarkers), '; Waiting for Image: ',filenamelist(counter+1,:)])    % plot a title onto the image

    end
end

save validx.dat validx -ascii -tabs
save validy.dat validy -ascii -tabs
msgboxwicon=msgbox('end.txt file detected, RTCorrCode stopped','Processing stopped!')
warning('end.txt file detected, RTCorrCode stopped')
end

%----------------------------------
%

function [processingtime]=fpstestfunc(grid_x,grid_y,filenamelist)
tic;

input_points_x=grid_x;
base_points_x=grid_x;

input_points_y=grid_y;
base_points_y=grid_y;

% [row,col]=size(base_points_x);      % this will determine the number of rasterpoints we have to run through
% [r,c]=size(filenamelist);                   % this will determine the number of images we have to loop through

base = uint8(mean(double(imread(filenamelist(1,:))),3));            % read in the base image ( which is always  image number one. You might want to change that to improve correlation results in case the light conditions are changing during the experiment
input = uint8(mean(double(imread(filenamelist(2,:))),3));       % read in the image which has to be correlated

input_points_for(:,1)=reshape(input_points_x,[],1);         % we reshape the input points to one row of values since this is the shape cpcorr will accept
input_points_for(:,2)=reshape(input_points_y,[],1);
base_points_for(:,1)=reshape(base_points_x,[],1);
base_points_for(:,2)=reshape(base_points_y,[],1);
input_correl(:,:)=cpcorr(input_points_for, base_points_for, input, base);           % here we go and give all the markers and images to process to cpcorr.m which ic a function provided by the matlab image processing toolbox
input_correl_x=input_correl(:,1);                                       % the results we get from cpcorr for the x-direction
input_correl_y=input_correl(:,2);                                       % the results we get from cpcorr for the y-direction

processingtime=toc;