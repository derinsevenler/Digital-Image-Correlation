% written by Chris

function [validx, validy]=pickpeak;

clear all

[Image,PathImage] = uigetfile('*.tif','Open Image');
cd(PathImage);
load('filenamelist');
filenumber=size(filenamelist);
filenumber=filenumber(1);

% filelist_generator

I=imread(Image); % read Image
Itemp=mean(double(I),3);
I=Itemp;
[Isizey, Isizex]=size(I);

prompt = {'How many peaks do you want to follow?'};
dlg_title = 'Manual peak picking';
num_lines= 1;
def     = {'2'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
numberofpeaks= str2num(cell2mat(answer(1,1)));

figure, imshow(uint8(I)); %show Image
axis('equal');
drawnow
title(sprintf('Mark the region of interest: Click on the on the lower left corner and and then on the upper right corner'));
hold on
cropxy(1,1:7)=0;

for i=1:numberofpeaks;
    
    [xprof, yprof]=ginput(2); % Get the Area of Interest
    cropxy(i,1) = i;
    cropxy(i,2) = (round(xprof(2,1)-xprof(1,1))/2)+xprof(1,1);
    cropxy(i,3) = (round(yprof(1,1)-yprof(2,1))/2)+yprof(2,1);
    cropxy(i,4) = xprof(1,1);
    cropxy(i,5) = yprof(2,1);
    cropxy(i,6) = round((xprof(2,1)-xprof(1,1))/2);
    cropxy(i,7) = round((yprof(1,1)-yprof(2,1))/2);
    % xmin = xprof(1,1);
    % xmax = xprof(2,1);
    % ymin = yprof(2,1);
    % ymax= yprof(1,1);
    
    plot(cropxy(i,2),cropxy(i,3),'o');
    drawnow;
    
end

tic; % start timer for time estimation
% I2 = imsubtract (I, imopen(I,strel('disk',15))); % subtract background
I2=I;
% [Isizey, Isizex]=size(I2);
% image(I2); %show with subtracted background5
% axis('equal');
t(1,1)=toc;
tic;

close all


% Start fitting process of the peaks, labeled by bwlabel

tic; % start timer
counter=0;
g = waitbar(0,'Processing image'); % nucleating the progress bar
fitcountertemp=size(cropxy); % number off peaks to cycle through
fitcounter=fitcountertemp(1,1); % number off peaks to cycle through
for c=1:fitcounter %start the loop to process all points
    waitbar(c/(fitcounter-1)); % growth of the progress bar
    cropI=imcrop(I2,[round(cropxy(c,4)) round(cropxy(c,5)) round(cropxy(c,6))*2 round(cropxy(c,7))*2]); % crop the region around the detected peak (bwlabel)
    
    % get line profile in x direction for the fitting routine
    
    xdata = [(round(cropxy(c,2))-round(cropxy(c,6))):1:(round(cropxy(c,2))+round(cropxy(c,6)))]; % x-coordinate for the fitting which is equivalent to the x coordinate in the image
    ydata=sum(cropI)/(2*cropxy(c,7)); % y-coordinate for the fitting which is equivalent to the greyvalues in the image - integrated in y direction of the image
    
    % fitting in x-direction
    % guess some parameters for the fitting routine --> bad guesses lead to
    % an error message which stops the fitting
    
    back_guess=(ydata(1)+ydata(round(cropxy(c,6))*2))/2; % guess for the background level - average of the first and last greyvalue
    sig1_guess=(cropxy(c,6)*2)/5; % guess for the peak width - take a fith of the cropping width
    amp_guess1=ydata(round((length(ydata))/2))-back_guess; % guess for the amplitude - take the greyvalue at the peak position
    mu1_guess=cropxy(c,2); % guess for the position of the peak - take the position from bwlabel
    
    % start fitting routine
    [x,resnormx,residual,exitflagx,output]  = lsqcurvefit(@gauss_onepk, [amp_guess1 mu1_guess sig1_guess back_guess], xdata, ydata); %give the initial guesses and data to the gauss fitting routine
    
    % show the fitting results
    xtest = [(round(cropxy(c,2))-round(cropxy(c,6))):0.1:(round(cropxy(c,2))+round(cropxy(c,6)))]; % x values for the plot of the fitting result
    ytest = (x(1)*exp((-(xtest-x(2)).^2)./(2.*x(3).^2))) + x(4); % y values of the fitting result
    yguess=(amp_guess1*exp((-(xtest-mu1_guess).^2)./(2.*sig1_guess.^2))) + back_guess; %y values for the guess plot
        plot(xdata,ydata,'o') % plot the experimental data
        hold on
        plot(xtest,ytest,'r') % plot the fitted function
        plot(xtest,yguess,'b') % plot the guessed function    
        drawnow
        hold off
%         pause
    % fitting in y-direction
    % guess parameters for the fitting routine --> bad guesses lead to
    % an error message which stops the fitting
    
    % get line profile in x direction for the fitting routine
    
    xdata = [(round(cropxy(c,3))-round(cropxy(c,7))):1:(round(cropxy(c,3))+round(cropxy(c,7)))]; % x data in y direction of the image
    ydata=sum(cropI')/(2*cropxy(c,6)); % integrate greyvalues in x direction and normalize it to the number of integrated lines
    
    % fitting in y-direction
    % guess parameters for the fitting routine --> bad guesses lead to
    % an error message which stops the fitting
    
    back_guess=(ydata(1)+ydata(round(cropxy(c,7))*2))/2; % guess for the background level - average of the first and last greyvalue
    sig1_guess=(cropxy(c,6)*2)/5; % guess for the peak width - take a fith of the cropping width
    amp_guess1=ydata(round((length(ydata))/2))-back_guess; % guess for the amplitude - take the greyvalue at the peak position
    mu1_guess=cropxy(c,3); % guess for the position of the peak - take the position from bwlabel
    
    % start fitting routine
    [y,resnormy,residual,exitflagy,output]  = lsqcurvefit(@gauss_onepk, [amp_guess1 mu1_guess sig1_guess back_guess], xdata, ydata); %give the initial guesses and data to the gauss fitting routine
    
    % show the fitting results
    xtest = [(round(cropxy(c,3))-round(cropxy(c,7))):0.1:(round(cropxy(c,3))+round(cropxy(c,7)))]; % x values for the plot of the fitting result
    ytest = (y(1)*exp((-(xtest-y(2)).^2)./(2.*y(3).^2))) + y(4); % y values of the fitting result
    yguess=(amp_guess1*exp((-(xtest-mu1_guess).^2)./(2.*sig1_guess.^2))) + back_guess; %y values for the guess plot
        plot(xdata,ydata,'o') % plot the experimental data
        hold on
        plot(xtest,ytest,'g') % plot the fitted function
        plot(xtest,yguess,'b') % plot the guessed function    
        drawnow
        hold off
%         pause
    % sort out the bad points and save the good ones in fitxy 
    % this matrix contains the to be used points from the first image
    
    if exitflagx>0 % if the fitting routine didn't find end before the 4000th iteration (check that in lsqcurvefit.m) then exitflag will be equal or smaller then 0
        cropxy(c,8)=1
        if exitflagy>0 % the same for the y direction fitting
            cropxy(c,8)=2
            if x(3)>0.05 % the width of the peak should be wider than 1 pixel - this is negotiable: different powder particle or cameras can give back results with very narrow peaks
                cropxy(c,8)=3
                if y(3)>0.05 % the same for y direction fitting
                    cropxy(c,8)=4
                    if resnormx<5000 % A measure of the "goodness" of fit is the residual, the difference between the observed and predicted data. (in the help file: Mathematics: Data Analysis and Statistics: Analyzing Residuals)
                        cropxy(c,8)=5
                        if resnormy<5000 % the same for the y- direction - - - a good value is as far as I know until now between 30 and 50. The good fits stay well beyond that (between 0 and 10)
                            cropxy(c,8)=6
                            if (round(x(2))-round(cropxy(c,6)))>0
                                cropxy(c,8)=7
                                if (round(x(2))+round(cropxy(c,6)))<Isizex
                                    cropxy(c,8)=8
                                    if (round(y(2))-round(cropxy(c,7)))>0
                                        cropxy(c,8)=9
                                        if (round(y(2))+round(cropxy(c,7)))<Isizey
                                            cropxy(c,8)=10
                                            counter=counter+1; 
                                            fitxy(counter,1)=c; % points  get their final number 
                                            fitxy(counter,2)=x(1); % fitted amplitude x-direction
                                            fitxy(counter,3)=abs(x(2)); % fitted position of the peak x-direction
                                            fitxy(counter,4)=abs(x(3)); % fitted peak width in x-direction
                                            fitxy(counter,5)=(x(4)); % fitted background in x-direction
                                            fitxy(counter,6)=y(1); % fitted amplitude y-direction
                                            fitxy(counter,7)=abs(y(2)); % fitted position of the peak y-direction
                                            fitxy(counter,8)=abs(y(3)); % fitted peak width in y-direction
                                            fitxy(counter,9)=abs(y(4)); % fitted background in y-direction
                                            fitxy(counter,10)=cropxy(c,6); % cropping width in x-direction
                                            fitxy(counter,11)=cropxy(c,7); % cropping width in ydirection
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end


close(g) % close progress bar window
t(1,3)=toc; % stop timer
image_time_s=t(1,3); % take time per image
estimated_totaltime_h=image_time_s*filenumber/3600; % calculate estimated time
sum(t);
total_time_h=sum(t);
close all

% plot image with peaks labeled by bwlabel (crosses) and the chosen points
% which are easy to fit with a gaussian distribution (circles)

figure, image(I2); %show Image
title(['Number of selected Images: ', num2str(filenumber), '; Estimated time [h] ', num2str((round(estimated_totaltime_h*10)/10)), ' Crosses are determined peaks, circles are chosen for  the analysis. If you want to run the analysis hit ENTER'])
axis('equal');
hold on;
plot(cropxy(:,2),cropxy(:,3),'+','Color','white') % peaks from bwlabel
plot(fitxy(:,3),fitxy(:,7),'o','Color','white'); % "good" points
drawnow

total_progress=1/filenumber;

pause

close all
fitlength=size(fitxy);
fitcounter=fitlength(1,1)
% again for all images
for m=1:(filenumber-1) % loop through all images
    tic; %start timer
    counter=0;
    f = waitbar(0,'Working on Image');
    I = imread(filenamelist(m,:)); %read image
    Itemp=mean(double(I),3);
    I=Itemp;
    
    
    % loop number
    for c=1:fitcounter %loop trough all points
        waitbar(c/(fitcounter-1)); %progress bar
        
        % load variables
        pointnumber=fitxy(c,(m-1)*12+1);
        amp_guess_x=fitxy(c,(m-1)*12+2);
        mu_guess_x=fitxy(c,(m-1)*12+3);
        sig_guess_x=fitxy(c,(m-1)*12+4);
        back_guess_x=fitxy(c,(m-1)*12+5);
        amp_guess_y=fitxy(c,(m-1)*12+6);
        mu_guess_y=fitxy(c,(m-1)*12+7);
        sig_guess_y=fitxy(c,(m-1)*12+8);
        back_guess_y=fitxy(c,(m-1)*12+9);
        crop_x=fitxy(c,(m-1)*12+10);
        crop_y=fitxy(c,(m-1)*12+11);
        
        % crop the area around the point to fit
        
        cropedI=imcrop(I,[(round(mu_guess_x)-round(crop_x)) (round(mu_guess_y)-round(crop_y)) 2*round(crop_x) 2*round(crop_y)]);
%         cropI=imsubtract (cropedI, imopen(cropedI,strel('disk',15))); % subtract background
        cropI=cropedI;%         imshow(cropI)
        % get line profile in x direction
        xdatax = [(round(mu_guess_x)-round(crop_x)):1:(round(mu_guess_x)+round(crop_x))];
        ydatax=sum(cropI)/(2*(crop_y));
        xguessx = [(round(mu_guess_x)-round(crop_x)):0.1:(round(mu_guess_x)+round(crop_x))];
        yguessx = (amp_guess_x*exp((-(xguessx-mu_guess_x).^2)./(2.*sig_guess_x.^2))) + back_guess_x;
        [x,resnormx,residualx,exitflagx,output]  = lsqcurvefit(@gauss_onepk, [amp_guess_x mu_guess_x sig_guess_x back_guess_x], xdatax, ydatax);
        xtestx = [(round(mu_guess_x)-round(crop_x)):0.1:(round(mu_guess_x)+round(crop_x))];
        ytestx = (x(1)*exp((-(xtestx-x(2)).^2)./(2.*x(3).^2))) + x(4);
%         plot(xdatax,ydatax,'o')
%         hold on
%         plot(xtestx,ytestx,'r')
%         plot(xguessx,yguessx,'b')
%         title(['Filename: ',filenamelist(m,:), '; Progress [%]: ',num2str((round(total_progress*10))/10), '; Tot. t [h] ', num2str((round(total_time_h*10)/10)), '; Est. t [h] ', num2str((round(estimated_totaltime_h*10)/10))])
%         drawnow
%         hold off
%         pause
        xdatay = [(round(mu_guess_y)-round(crop_y)):1:(round(mu_guess_y)+round(crop_y))];
        ydatay=sum(cropI')/(2*(crop_y));
        xguessy = [(round(mu_guess_y)-round(crop_y)):0.1:(round(mu_guess_y)+round(crop_y))];
        yguessy = (amp_guess_y*exp((-(xguessy-mu_guess_y).^2)./(2.*sig_guess_y.^2))) + back_guess_y;
        [y,resnormy,residualy,exitflagy,output]  = lsqcurvefit(@gauss_onepk, [amp_guess_y mu_guess_y sig_guess_y back_guess_y], xdatay, ydatay);
        xtesty = [(round(mu_guess_y)-round(crop_y)):0.1:(round(mu_guess_y)+round(crop_y))];
        ytesty= (y(1)*exp((-(xtesty-y(2)).^2)./(2.*y(3).^2))) + y(4);
%         plot(xdatay,ydatay,'o')
%         hold on
%         plot(xtesty,ytesty,'g')
%         plot(xguessy,yguessy,'b')
%         title(['Filename: ',filenamelist(m,:), '; Progress [%]: ',num2str((round(total_progress*10))/10), '; Tot. t [h] ', num2str((round(total_time_h*10)/10)), '; Est. t [h] ', num2str((round(estimated_totaltime_h*10)/10))])
%         drawnow
%         hold off
%         pause
        
        if exitflagx>0
            if exitflagy>0
                counter=counter+1;
                fitxy(counter,m*12+1)=pointnumber;
                fitxy(counter,m*12+2)=abs(x(1));
                fitxy(counter,m*12+3)=abs(x(2));
                fitxy(counter,m*12+4)=abs(x(3));
                fitxy(counter,m*12+5)=abs(x(4));
                fitxy(counter,m*12+6)=abs(y(1));
                fitxy(counter,m*12+7)=abs(y(2));
                fitxy(counter,m*12+8)=abs(y(3));
                fitxy(counter,m*12+9)=abs(y(4));
                fitxy(counter,m*12+10)=crop_x;
                fitxy(counter,m*12+11)=crop_y;
                fitxy(counter,m*12+12)=resnormx;
                
            end
        end
        
        
    end
    imshow(uint8((I))); %show Image
    title(['Filename: ',filenamelist(m,:), '; Progress [%]: ',num2str((round(total_progress*10))/10), '; Tot. t [h] ', num2str((round(total_time_h*10)/10)), '; Est. t [h] ', num2str((round(estimated_totaltime_h*10)/10))])    
    axis('equal');
    hold on;
    plot(fitxy(:,m*12+3),fitxy(:,m*12+7),'o','Color','white'); % "good" points
    drawnow
    
    total_progress=1/filenumber;
    
    %         pause
    hold off;
    
%     plot(fitxy(:,m*12+1),fitxy(:,m*12+12),'+');
%     title(['Filename: ',filenamelist(m,:), '; Progress [%]: ',num2str((round(total_progress*10))/10), '; Tot. t [h] ', num2str((round(total_time_h*10)/10)), '; Est. t [h] ', num2str((round(estimated_totaltime_h*10)/10))])
    fitcounter=counter;
    close(f);
    time(m)=toc;
    total_time_s=sum(time);
    total_time_h=sum(time)/3600;
    image_time_s=total_time_s/m;
    estimated_totaltime_h=image_time_s*(filenumber)/3600;
    progress_percent=total_time_h/estimated_totaltime_h*100;
    total_progress=(m+1)/(filenumber)*100;
    
end  

% save the stuff
save fitxy.dat fitxy -ascii -tabs

[validx,validy]=sortvalidpoints(fitxy);
title(['Processing Images finished!'])

save('validx');
save('validy');

save validx.dat validx -ascii -tabs;
save validy.dat validy -ascii -tabs;

close all
