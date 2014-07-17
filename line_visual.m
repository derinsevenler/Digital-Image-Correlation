function line_visual(filenamelist,raw_peak_results,int_width_y_mark);

%"Stand alone" code to visualize DIC .dat results files
%Programmed by Rob slightly changed by Chris
%Last revised: 2//06


%Load the necessary files

if exist('filenamelist')==0
load('filenamelist')
end
if exist('raw_peak_results')==0
load('raw_peak_results.dat')
end
if exist('int_width_y_mark')==0
load('int_width_y_mark.dat')
end

%Number of images to analyze is r; also initialize xpos, y_mark

[r,c]=size(filenamelist);
xpos_o=raw_peak_results(1,:);
y_mark=int_width_y_mark(2,1);

%Start analysis of results

g=waitbar(0,'Processing the results...');

for m=2:(r-1)
    
    %Calculate displacements and engineering strains (i.e. average strains between
    %the two selected points...similar to an ISDG measurement)
    
    waitbar(m/(r-1));
    warning off MATLAB:divideByZero
    
    xdis=raw_peak_results(m,:)-xpos_o;  
    epsxx=(xdis(2)-xdis(1))/(xpos_o(2)-xpos_o(1));     
    
    epsxx_plot(1)=0;
    epsxx_plot(m+1)=epsxx*100;
    
    %Initialize and then save displacement and strain field data for each image
    
    xdis_all(m,:)=xdis;      
    epsxx_all(m)=epsxx;
end

epsxx_all=epsxx_all';

close(g)

save disp_x.dat xdis_all -ascii -tabs 
save strain_x.dat epsxx_all -ascii -tabs

%Show user strain vs. image number plot and save the data and figure

l=figure;
plot(1:r,epsxx_plot,'.-b')
axis([1 (max(r)*1.1) 0 (max(epsxx_plot)*1.1)])
xlabel('Image Number (Image #1 defined to be the very first image loaded for analysis)')
ylabel('Engineering Strain [%]')
title('Engineering Strain vs. Image Number') 

strain_image_data_x=[(1:r)' epsxx_plot'];
save strain_image_x.txt strain_image_data_x -ascii -tabs
saveas(l,'strain_image_x.fig')

%Ask user if a full-figure stress vs. strain plot is desired; save data and figure too

plot_selection = menu(sprintf('The displacement and strain values from each image have been saved and plotted.\n  An Eng Stress vs. Eng Strain curve can now be generated if the files "time_stress.txt" (raw data from Labview, stress values in [MPa]) and "time_image.txt" (saved by AutoIt)\n are available in the image directory.  Matlab will sort through these files and select the correct stress value for a given image and calculated strain using system time as the reference.\n  (Note: The matching logic used in this step to correlate the stress and strain values requires that the raw data aquisition rate be at least 2 data points per second.)\n  Click "Generate Plot" to do this or click "Exit" to quit Matlab.'),'Generate Plot','Exit');

if plot_selection == 1
    
    stress_strain_match
    
    n=figure;
    plot(epsxx_plot,stress_image(:,2),'.-b')
    axis([0 (max(epsxx_plot)*1.1) 0 (max(stress_image(:,2))*1.1)])
    xlabel('Engineering Strain [%]')
    ylabel('Engineering Stress [MPa]')
    title('Engineering Stress vs. DIC Engineering Strain') 
    
    stress_strain_data_x=[epsxx_plot' stress_image(:,2)];
    save stress_strain_x.txt stress_strain_data_x -ascii -tabs
    saveas(n,'stress_strain_x.fig')
    
elseif plot_selection == 2
    return
end

%Ask user if a figure with two subplots is desired...good for presentation animations

plot_selection2 = menu(sprintf('A figure with two subplots can now be generated for use in animated presentations.\n  (Caution: The .avi movie saved in this step can become very large depending on the number of images analyzed.)\n  Click "Generate Subplots" to do this or click "Exit" to quit Matlab.'),'Generate Subplots','Exit');

if plot_selection2 == 1
    
    %Define counter for image capture and load stress data
    
    mov_count=0;
    
    for M=1:r
        
        subplot(2,1,1)
        imshow(filenamelist(M,:))
        hold on 
        plot(raw_peak_results(1,1),y_mark,'xr','markersize',10)           %,'markersize',12
        plot(raw_peak_results(1,2),y_mark,'xr','markersize',10)
        
        if M > 1
            
            plot(raw_peak_results(M,1),y_mark,'xg','markersize',10)
            plot(raw_peak_results(M,2),y_mark,'xg','markersize',10)
        end
        
        title(['Image Correlation Results',sprintf(' (Current image: %10s)',filenamelist(M,:))],'fontweight','bold')
        drawnow
        hold off
        
        subplot(2,1,2)
        hold on
        plot(epsxx_plot(M),mov_stress(M),'og')
        box on
        axis([0 (max(epsxx_plot)*1.1) 0 (max(mov_stress)*1.1)])
        xlabel('Engineering Strain [%]','fontsize',8,'fontweight','bold')
        ylabel('Engineering Stress [MPa]','fontsize',8,'fontweight','bold')
        title('Engineering Stress vs. DIC Engineering Strain','fontweight','bold')
        drawnow
        hold off
        
        
    end
    
elseif plot_selection2 == 2
    return
end



