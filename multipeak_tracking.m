%% written by Chris

% updated last 9/8/2010
% this function allows to track multiple peaks along one axis. The
% resolution can be much higher than the cross correlation methods of
% digital image correlation.

% function goes here
function [validx,validy]=multipeak_tracking;

%% Load image and ask for orientation and pos/neg

if exist('filenamelist')==0
    load('filenamelist')            % file with the list of filenames to be processed
end

filenumber=size(filenamelist);
filenumber=filenumber(1);

h=figure;
imshow(filenamelist(1,:))           % show the first image
title('First Image')        % put a title

ImageA=(uint8(mean(double(imread(filenamelist(1,:))),3)));

sel_pos_neg = menu(sprintf('Do you want use the negative image for analysis?'),...
    'As is', 'Negativ', 'Cancel')
if sel_pos_neg==3
    close all;
    return
end
if sel_pos_neg==2
    ImageA=255-ImageA;
    imshow(ImageA)           % show the first image
    negativ=1;
end
if sel_pos_neg==1
    negativ=0;
end

sel_horz_vert = menu(sprintf('Do you want to measure horizontal or vertical?'),...
    'Horizontal', 'Vertical', 'Cancel')
if sel_horz_vert==3
    close all;
    return
end
if sel_horz_vert==2
    ImageA=ImageA';
    imshow(ImageA)           % show the first image
    orientation=90;
end
if sel_horz_vert==1
    orientation=0;
end

[imr,imc]=size(ImageA);


%% select path for displacement extraction along markers
%Select point in middle of marker
xlabel('Location On Image [Pixels]')
ylabel('Location On Image [Pixels]')
title(['Click on the center of the sample. Width: ', num2str(imc),...
    '; Height: ',num2str(imr)])
marker_pt=round(ginput(1));
x_mark = marker_pt(1);
y_mark = marker_pt(2);
hold on
line([1 imc],[y_mark y_mark],'Color','r')


%Prompt for integration width
prompt = {'Enter integration width [Pixels]:'};
dlg_title = 'Input integration width for the analysis';
num_lines= 1;
def     = {'20'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
int_width = str2num(cell2mat(answer(1)));

line([1 imc],[y_mark-int_width/2 y_mark-int_width/2],'Color','g')
line([1 imc],[y_mark+int_width/2 y_mark+int_width/2],'Color','g')


%% Calculate line profile data and select peaks
xdata = [1:1:imc];
ydata= sum(ImageA((y_mark-int_width/2):(y_mark+int_width/2),:),1)/int_width;

g=figure;
plot(xdata,ydata)
xlabel('Location On Image [Pixels]')
ylabel('Pixel Intensity Value')

sel_peakselect = menu(sprintf('Do you want to select single or all peaks?')...
    ,'Single Select', 'All Select', 'Cancel')
if sel_peakselect==3
    close all;
    return
end

% If single peaks should be chosen
if sel_peakselect==1
    sel_peakselect_2=1
    i=1
    hold on
    while sel_peakselect_2==1
        title('Click on the left of the chosen peak.')
        marker_pt=round(ginput(1));
        x_peak_left = marker_pt(1);
        y_peak_left = marker_pt(2);
        line([x_peak_left x_peak_left], [1 max(ydata)],'Color','r')
        
        title('Click on the right of the chosen peak.')
        marker_pt=round(ginput(1));
        x_peak_right = marker_pt(1);
        y_peak_right = marker_pt(2);
        line([x_peak_right x_peak_right], [1 max(ydata)],'Color','r')
        
        plot(xdata(x_peak_left:x_peak_right),ydata(x_peak_left:x_peak_right),'k')
        peak_pos(i,:)=[x_peak_left x_peak_right];
        
        i=i+1;
        sel_peakselect_2 = menu(sprintf('Do you want to select another peak?')...
            ,'One more', 'No, lets process')
    end
end

% If all peaks should be found automatically
if sel_peakselect==2
    
end

%% Start first fitting
fitcountertemp=size(peak_pos); % number off peaks to cycle through
fitcounter=fitcountertemp(1,1); % number off peaks to cycle through

for c=1:fitcounter %start the loop to process all points
    fitdatax=xdata(peak_pos(c,1):peak_pos(c,2));
    fitdatay=ydata(peak_pos(c,1):peak_pos(c,2));
    
    % fitting in x-direction
    % guess some parameters for the fitting routine --> bad guesses lead to
    % an error message which stops the fitting
    
    back_guess=(ydata(peak_pos(c,2))+ydata(peak_pos(c,2)))/2; % guess for the background level - average of the first and last greyvalue
    sig1_guess=(peak_pos(c,2)-peak_pos(c,1))/2; % guess for the peak width - take half of the cropping width
    amp_guess1=max(fitdatay); % guess for the amplitude - take the greyvalue at the peak position
    mu1_guess=(peak_pos(c,2)+peak_pos(c,1))/2; % guess for the position of the peak - take the position from bwlabel
    
    % start fitting routine
    [x,resnormx,residual,exitflagx,output]  = lsqcurvefit(@gauss_onepk, [amp_guess1 mu1_guess sig1_guess back_guess], fitdatax, fitdatay); %give the initial guesses and data to the gauss fitting routine
    
    % show the fitting results
    xtest = fitdatax; % x values for the plot of the fitting result
    ytest = (x(1)*exp((-(xtest-x(2)).^2)./(2.*x(3).^2))) + x(4); % y values of the fitting result
    yguess=(amp_guess1*exp((-(xtest-mu1_guess).^2)./(2.*sig1_guess.^2))) + back_guess; %y values for the guess plot
    plot(xtest,ytest,'r') % plot the fitted function
    plot(xtest,yguess,'b') % plot the guessed function
    drawnow
    
    peaks(c,1,:)=x;
end

sel_go = menu(sprintf('Do you want proceed?')...
    ,'Do it!', 'No, stop!')

if sel_peakselect==2
    close all;
    return
end

cutout=peak_pos(:,2)-peak_pos(:,1)

for m=1:(filenumber-1) % loop through all images
    ImageA=((mean(double(imread(filenamelist(m,:))),3)));
    
    if sel_horz_vert==2
        ImageA=ImageA';
        imshow(ImageA)           % show the first image
        orientation=90;
    end
    if sel_pos_neg==2
        ImageA=255-ImageA;
        imshow(ImageA)           % show the first image
        negativ=1;
    end
    
    image=m
    
    xdata = [1:1:imc];
    ydata= sum(ImageA((y_mark-int_width/2):(y_mark+int_width/2),:),1)/int_width;
    
    plot(xdata,ydata)
    title(['Image Number: ', num2str(m)])
    hold on
    
    
    for c=1:fitcounter %start the loop to process all points
        
        if (peaks(c,m,2)+cutout(c)/2)>imc
        fitdatax=xdata((imc-cutout(c)):imc);
        fitdatay=ydata((imc-cutout(c)):imc);
        else
        fitdatax=xdata((peaks(c,m,2)-cutout(c)/2):(peaks(c,m,2)+cutout(c)/2));
        fitdatay=ydata((peaks(c,m,2)-cutout(c)/2):(peaks(c,m,2)+cutout(c)/2));
        end
        
        % fitting in x-direction
        % guess some parameters for the fitting routine --> bad guesses lead to
        % an error message which stops the fitting
        
        back_guess=peaks(c,m,4); % guess for the background level - average of the first and last greyvalue
        sig1_guess=peaks(c,m,3); % guess for the peak width - take half of the cropping width
        amp_guess1=peaks(c,m,1); % guess for the amplitude - take the greyvalue at the peak position
        mu1_guess=peaks(c,m,2); % guess for the position of the peak - take the position from bwlabel
        
        % start fitting routine
        [x,resnormx,residual,exitflagx,output]  = lsqcurvefit(@gauss_onepk, [amp_guess1 mu1_guess sig1_guess back_guess], fitdatax, fitdatay); %give the initial guesses and data to the gauss fitting routine
        % show the fitting results
        xtest = fitdatax; % x values for the plot of the fitting result
        ytest = (x(1)*exp((-(xtest-x(2)).^2)./(2.*x(3).^2))) + x(4); % y values of the fitting result
        yguess=(amp_guess1*exp((-(xtest-mu1_guess).^2)./(2.*sig1_guess.^2))) + back_guess; %y values for the guess plot
        plot(xtest,ytest,'r') % plot the fitted function
        plot(xtest,yguess,'g') % plot the guessed function
        drawnow
        peaks(c,m+1,:)=x;
    end
    hold off
    Amplitude=peaks(:,m,1)';
    Peakposition=peaks(:,m,2)';
    Peakwidth=peaks(:,m,3)';
    Background=peaks(:,m,4)';
    dlmwrite('backupPeakposition.txt', Peakposition' , 'delimiter', '\t', '-append');       % Here we save the result from each image; if you are desperately want to run this function with e.g. matlab 6.5 then you should comment this line out. If you do that the data will be saved at the end of the correlation step - good luck ;-)
    dlmwrite('backupAmplitude.txt', Amplitude' , 'delimiter', '\t', '-append');       % Here we save the result from each image; if you are desperately want to run this function with e.g. matlab 6.5 then you should comment this line out. If you do that the data will be saved at the end of the correlation step - good luck ;-)
    dlmwrite('backupPeakwidth.txt', Peakwidth' , 'delimiter', '\t', '-append');       % Here we save the result from each image; if you are desperately want to run this function with e.g. matlab 6.5 then you should comment this line out. If you do that the data will be saved at the end of the correlation step - good luck ;-)
    dlmwrite('backupBackground.txt', Background' , 'delimiter', '\t', '-append');       % Here we save the result from each image; if you are desperately want to run this function with e.g. matlab 6.5 then you should comment this line out. If you do that the data will be saved at the end of the correlation step - good luck ;-)
end

%% Save data
validx=peaks(:,:,2);
validy=peaks(:,:,3);

Amplitude=peaks(:,:,1)';
Peakposition=peaks(:,:,2)';
Peakwidth=peaks(:,:,3)';
Background=peaks(:,:,4)';

save peaks;
save Amplitude.dat Amplitude -ascii -tabs;
save Peakposition.dat Peakposition -ascii -tabs;
save Peakwidth.dat Peakwidth -ascii -tabs;
save Background.dat Background -ascii -tabs;
save validx_mp.dat validx -ascii -tabs;
save validy_mp.dat validy -ascii -tabs;