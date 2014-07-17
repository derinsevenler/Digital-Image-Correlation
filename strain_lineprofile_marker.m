function [validx,validy]=strain_lineprofile_marker;

%Code to use ISDG markers for 1-D Matlab DIC strain measurement
%Programmed by Dan and Rob
%Changed by Chris
%Last revision: 3. February 2008


%Load the filenames to analyze
load('filenamelist')
[r,c]=size(filenamelist);

%Open base image
[FileNameBase,PathNameBase] = uigetfile( ...
    {'*.bmp;*.tif;*.jpg;*.TIF;*.BMP;*.JPG','Image files (*.bmp,*.tif,*.jpg)';'*.*',  'All Files (*.*)'}, ...
    'Open base image for marker selection (color image best)');
cd(PathNameBase);
im_base = imread(FileNameBase);
im_basetemp=mean(double(im_base),3);
im_base=im_basetemp;

maxcolor=max(max(im_base));
mincolor=min(min(im_base));
caxis([mincolor maxcolor]);
imagesc(im_base);

[imr,imc]=size(im_base);

%Select point in middle of marker
xlabel('Location On Image [Pixels]')
ylabel('Location On Image [Pixels]')
title(sprintf('Click on the center of one marker.'))
marker_pt=round(ginput(1));
x_mark = marker_pt(1);
y_mark = marker_pt(2);

%Prompt for integration width
prompt = {'Enter integration width [Pixels]:'};
dlg_title = 'Input integration width for the analysis';
num_lines= 1;
def     = {'40'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
int_width = str2num(cell2mat(answer(1)));

%Calculate line profile data
xdata = [1:1:imc];
ydata= sum(im_base((y_mark-int_width/2):(y_mark+int_width/2),:),1)/int_width;

h=figure;
set(h,'Position',[150,50,800,600])
plot(xdata,ydata)
xlabel('Location On Image [Pixels]')
ylabel('Pixel Intensity Value')
title(sprintf('First, single click at the base and middle of the LEFT PEAK (vertically even with the flat data).\n  Then, single click on the max point of the RIGHT PEAK.'))
[xprof, yprof]=ginput(2);
mu1_guess = xprof(1,1);
mu2_guess = xprof(2,1);
back_guess = yprof(1,1);
amp_guess1 = yprof(2,1);
amp_guess2 = yprof(2,1);
sig1_guess = 3;
sig2_guess = 3;

%Save figure
saveas(h,'center_select.fig')

%Close all open images
close all

%Initialize status bar and figure
k=figure;
g=waitbar(0,'Processing the images...');
set(g,'Position',[275,50,275,50])  %[left bottom width height]
% results_prof = zeros(r,2);
tic;
% peaklength=40;
validx=zeros(2,r);
validy=zeros(2,r);

for i=1:(r-1)
    
    waitbar(i/(r-1));
    
    %Read in images from filenamelist
    im_input = mean(double(imread(filenamelist(i,:))),3);

    %Generate greyscale data (y-values)
    ydata= sum(im_input((int8(y_mark-int_width/2)):(int8(y_mark+int_width/2)),:),1)/int_width;
    warning off all
    
    %Fit gaussian function to peaks
    [x,resnorm,residual,exitflag,output]  = lsqcurvefit(@gauss_twopk, [amp_guess1 mu1_guess sig1_guess amp_guess2 mu2_guess sig2_guess back_guess], xdata, ydata);
    
    %plot fit and data
    xtest = [0:0.1:imc];
    ytest = (x(1)*exp((-(xtest-x(2)).^2)./(2.*x(3).^2))) + (x(4)*exp((-(xtest-x(5)).^2)./(2.*x(6).^2))) + x(7);
    plot(xdata,ydata,'.')
    axis([0 imc 0 (max(ydata)*1.1)])
    xlabel('Location On Image [Pixels]')
    ylabel('Pixel Intensity Value')
    title('Pixel Intensity Value Fit vs. Location On Image') 
    hold on;
    plot(xtest,ytest,'r')
    hold off;
    
    %Write results to matrix
    results_prof(i+1,1) = x(2);
    results_prof(i+1,2) = x(5);
    results_prof(1,1) = results_prof(2,1);
    results_prof(1,2) = results_prof(2,2);
    
    validx(:,i)=[x(2) x(5)];
    validy(:,i)=[y_mark y_mark];

    %Update guesses
    amp_guess1 = x(1); mu1_guess = x(2); sig1_guess = x(3); amp_guess2 = x(4); mu2_guess = x(5); sig2_guess = x(6); back_guess = x(7);
    
end
toc
close(g)

%Save figure
saveas(k,'peak_fit.fig')

%Save locations of two peaks for all images and save integration width
save raw_peak_results.dat results_prof -ascii -tabs
save int_width_y_mark.dat int_width y_mark -ascii -tabs
save validx_line_profile.dat validx -ascii -tabs
save validy_line_profile.dat validy -ascii -tabs

line_visual(filenamelist,results_prof,[int_width; y_mark]);