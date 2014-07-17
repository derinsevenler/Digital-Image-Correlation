
% Initialize data
% written by Chris and Dan

% Displacement.m allows you to analyze the data you aquiered with the
% correlation, fitting or mean routine. It only needs the validx and
% validy and can calculate strain from it. Before you start you should 
% consider cleaning up the data as described in the guide. After that step
% you can analyze parts of your data, or the full set. Try to use also the
% console command, e.g. if you want to analyze only image 100-110 since
% something really interesting happend there, load validx and validy into
% your workspace and call
% displacement(validx(:,100:110),validy(:,100:110));
% In this case displacement only loads the important images and you can
% clean this part of your data set.

% Changed 3. February 2008


function [validx,validy]=displacement(validx,validy);

%load data in case you did not load it into workspace yet
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

%define the size of the data set
sizevalidx=size(validx);
sizevalidy=size(validy);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);

%calculate the displacement relative to the first image in x and y
%direction
clear displx;
validxfirst=zeros(size(validx));
validxfirst=mean(validx(:,1),2)*ones(1,sizevalidx(1,2));
displx=validx-validxfirst;
clear validxfirst
clear disply;
validyfirst=zeros(size(validy));
validyfirst=mean(validy(:,1),2)*ones(1,sizevalidy(1,2));
disply=validy-validyfirst;
clear validyfirst

%Prompt user for type of plotting / visualization
selection10 = menu(sprintf('How do you want to visualize your data?'),'3D Mesh Plot of Displacement'...
    ,'Full Strain Plots','Strain Measurement between 2 Points','1D Average Strain Measurement',...
    'Rotate Orientation (exchange x and y)','Remove badly tracked marker, one by one (Position)',...
    'Delete multible markers (Position)','Delete markers from displacement vs. position plot',...
    'Delete points moving relative to their neighbours','Select Markers to Analyze ',...
    'Save validx and validy','Average a couple of images','Cancel');

% Selection for Cancel - All windows will be closed and you jump back to
% the command line
if selection10==13
    close all;
    return
end

% This selection will average up a specified number of images to reduce the
% noise of the data set. I would like to point out that you will need to
% average your other sensor data (e.g. load data), too, to match it to your
% strain data.
if selection10==12
    prompt = {'How many images would you like to combine as a base image?'};
    dlg_title = 'Input number of images:';
    num_lines= 1;
    def     = {'5'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if str2num(cell2mat(answer(1)))==0
        disp('Get out, you changed your mind?')
        [validx validy]=displacement(validx,validy);
        return
    else
        baseimages = str2num(cell2mat(answer(1)));
        if baseimages==[]
            disp('Give me a number, will you?')
            [validx validy]=displacement(validx,validy);
            return
        end
        if baseimages>loopimages
            disp('That is too large?!')
        else
            baseimagemean=mean(validx(:,1:baseimages),2);
            validx(:,1:baseimages-1)=[];
            validx(:,1)=baseimagemean;
            baseimagemean=mean(validy(:,1:baseimages),2);
            validy(:,1:baseimages-1)=[];
            validy(:,1)=baseimagemean;
        end
    end
    [validx validy]=displacement(validx,validy);
    return
end

% Save validx and validy, very useful if you cleaned up your dataset. Data
% will be saved as -ascii text file. If you send data like this by email
% you can reduce the size tremendously by compressing it. Use ZIP or RAR.
if selection10==11
    [FileName,PathName] = uiputfile('validx_corr.dat','Save validx');
    if FileName==0
        disp('You did not save your file!')
        [validx validy]=displacement(validx,validy);
        return
    else
        cd(PathName)
        save(FileName,'validx','-ascii')
        [FileName,PathName] = uiputfile('validy_corr.dat','Save validy');
        if FileName==0
            disp('You did not save your file!')
            [validx validy]=displacement(validx,validy);
        else
            cd(PathName)
            save(FileName,'validy','-ascii')
        end
        [validx validy]=displacement(validx,validy);
        return
    end
end

% Select Points from detailed Analysis
if selection10==10
    [validx validy validxbackup validybackup]=ppselection_func(validx,validy);
    if validx==0
        validx=validxbackup;
        validy=validybackup;
    end
    if validy==0
        validx=validxbackup;
        validy=validybackup;
    end
    [validx validy]=displacement(validx,validy);
end

% Remove markers moving relativ to their neighbours
if selection10==9
    [validx,validy,displx,disply]=delete_jumpers(validx,validy,displx,disply);
    [validx validy]=displacement(validx,validy);
end

% Remove markers from the displacement vs. position plot
if selection10==8
    [validx,validy,displx,disply]=removepoints_func(validx,validy,displx,disply);
    [validx validy]=displacement(validx,validy);
end

% Remove bad points
if selection10==7
    [validx,validy]=removepoints_func2(validx,validy);
    [validx validy]=displacement(validx,validy);
end

% Remove bad points
if selection10==6
    [validx validy]=removepoints_func3(validx,validy);
    [validx validy]=displacement(validx,validy);
end

% Rotate Matrix
if selection10==5
    [validx, validy]=rotatematrix(validx,validy);
    [validx validy]=displacement(validx,validy);
end

% 1D Strain plot using average strains for ELASTIC STRAIN only
if selection10==4
    [validx validy]=strain_1D_average_func(validx, validy,displx,disply);
    [validx validy]=displacement(validx,validy);
end

% 1D Strain plot
if selection10==3
    [validx, validy,displx,disply]=strain_1D_2Points_func(validx, validy,displx,disply);
    [validx validy]=displacement(validx,validy);
end

% Fast plotting, cropping needed for polynomial fit
if selection10==2
    [validx, validy,displx,disply]=polyfit3D(validx, validy,displx,disply);
    [validx validy]=displacement(validx,validy);
end

% 3D Mesh Plotting
if selection10==1
    if sizevalidx(1,1)>2
    [validx, validy,displx,disply]=meshplot(validx,validy,displx,disply);
    else
        disp('You need at least three markers to display the 3D-plot')
        msgbox('You need at least three markers to display the 3D-plot','3D-Plot','warn');
    end
    [validx validy]=displacement(validx,validy);
end

%---------------------------------

function [validx,validy,displx,disply]=delete_jumpers(validx,validy,displx,disply);

% written by Chris

% This is a filter which helps to find jumpy data points which are 
% oscillating or stop moving.
% The Filter starts by finding the next 10 datapoint neighbours 
% (num_neighbours), calculates their mean position and then plots the
% difference between each data point and its neighbours versus image
% number. If a data point is jumping around it will show up as a spike. But
% be careful, one bad one will also affect his neighbours, therefore its
% worthwhile to use this filter step by step.

% Changed 3. February 2008

num_neighbours=10;

doitonemoretime=1

while doitonemoretime==1
    % defining variables
    sizevalidx=size(validx);
    sizevalidy=size(validy);
    looppoints=sizevalidx(1,1);
    loopimages=sizevalidx(1,2);

    % clear the used ones
    clear validxtemp
    clear validytemp
    clear meandistancetemp
    clear sizevalidxtemp
    clear sizevalidytemp
    clear looppointstemp
    clear loopimagestemp
    clear max_distance
    clear min_distance
    clear dist_matrix
    clear dist_sort
    clear dist_index
    clear meandistance
    tic
    % calculate the distance to the next data points
%     dist_matrix=zeros(looppoints,looppoints);
    meandistance=zeros(sizevalidx);
    
    g=waitbar(0,'Processing the markers...');
    for i=1:looppoints
        waitbar(i/looppoints);
        dist_matrix=(((validx(:,1)-validx(i,1)).^2+(validy(:,1)-validy(i,1)).^2).^(0.5))';
        %   end

        % find the next neighbours by indexing the ones closest
        [dist_sort, dist_index]=sort(dist_matrix);

        % take the mean position of the closest data points of each for all
        % images
        meandistance(i,:)= validx(i,:)-mean(validx(dist_index(2:num_neighbours),:),1);
        max_distance(i)= max(diff(meandistance(i,:)-meandistance(i,1)));
        min_distance(i)= min(diff(meandistance(i,:)-meandistance(i,1)));
    end
    close(g)
toc
    for i=1:looppoints
        plot(diff(meandistance(i,:)-meandistance(i,1)))
        hold on
    end
   toc     
    % Select an upper and lower boundary
    xlabel('Image number[ ]')
    ylabel('Relative marker dispacement [Pixels]')
    title(sprintf('Define the upper and lower bound by clicking above and below the valid points'))
    marker_pt=(ginput(1));
    x_mark(1,1) = marker_pt(1);
    y_mark(1,1) = marker_pt(2);
    plot([1;loopimages],[y_mark(1,1);y_mark(1,1)],'r');

    title(sprintf('Define the upper and lower bound by clicking above and below the valid points'))
    marker_pt=(ginput(1));
    x_mark(1,2) = marker_pt(1);
    y_mark(1,2) = marker_pt(2);
    plot([1;loopimages],[y_mark(1,2);y_mark(1,2)],'r');

    upperbound=max(y_mark);
    lowerbound=min(y_mark);

    hold off
    
    validxtemp=validx;
    validytemp=validy;
    meandistancetemp=meandistance;
    
    
    validxtemp(find(max_distance>upperbound | min_distance<lowerbound),:)=[];
    validytemp(find(max_distance>upperbound | min_distance<lowerbound),:)=[];
    meandistancetemp(find(max_distance>upperbound |min_distance<lowerbound),:)=[];
    sizevalidxtemp=size(validxtemp);
    sizevalidytemp=size(validytemp);
    looppointstemp=sizevalidxtemp(1,1);
    loopimagestemp=sizevalidxtemp(1,2);
    
    for i=1:looppointstemp
        plot(diff(meandistancetemp(i,:)-meandistancetemp(i,1)))
        hold on
    end
    plot([1;loopimages],[y_mark(1,1);y_mark(1,1)],'r');
    plot([1;loopimages],[y_mark(1,2);y_mark(1,2)],'r');
    
    hold off

    selection_filter = menu('Do you like the result?','Take it as is','Want to select more','Try again','Cancel');
    if selection_filter==1
        validx=validxtemp;
        validy=validytemp;
        doitonemoretime=0;
    elseif selection_filter==2
        validx=validxtemp;
        validy=validytemp;
        doitonemoretime=1;
    elseif selection_filter==3
        doitonemoretime=1;
    elseif selection_filter==4
        return
    end
end

%---------------------------------
% Rotate Matrix
% written by Chris
function [validx, validy]=rotatematrix(validx,validy);
validxrot=validx;
clear validx;
validyrot=validy;
clear validy;
validy=validxrot;
validx=validyrot;


%---------------------------------
% Delete points from the displacement plot
% written by Chris
function [validx,validy,displx,disply] = removepoints_func(validx,validy,displx,disply) ; %delete points

close all

if exist('validx')==0
    [validx,Pathvalidx] = uigetfile('*.mat; *.txt','Open validx.mat or validx.txt');
    cd(Pathvalidx);
    validx=importdata(validx,'\t');
    [validy,Pathvalidy] = uigetfile('*.mat;*.txt','Open validy.mat or validy.txt');
    cd(Pathvalidy);
    validy=importdata(validy,'\t');
end

selectremove1 = menu(sprintf('Do you want to delete makers?'),'Yes','No');
if selectremove1==2
    return
end

% if yes
if selectremove1==1

    sizevalidx=size(validx);
    sizevalidy=size(validy);

    selectionremove2=selectremove1;
    counter=0
    sizevalidx=size(validx);
    looppoints=sizevalidx(1,1);
    loopimages=sizevalidx(1,2);
    defaultimage=loopimages
    numberbadpoints=0

    while selectionremove2==1
        counter=counter+1
        clear xplot
        clear sizevalidx
        clear selectremove11
        clear selection2
        %         clear badpoints

        sizevalidx=size(validx);
        looppoints=sizevalidx(1,1);
        loopimages=sizevalidx(1,2);

        % update temporary matrices
        validxtemp=validx;
        validytemp=validy;
        displxtemp=displx;
        displytemp=disply;

        % get the image number from which the bad points will be chosen
        prompt = {'From which image do you want to delete markers?'};
        dlg_title = 'Marker removal';
        num_lines= 1;
        if numberbadpoints==0
            defaultimage=loopimages
        end
        if numberbadpoints~0
            defaultimage=numberbadpoints
        end
        def     = {num2str(defaultimage)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        numberbadpoints = str2num(cell2mat(answer(1,1)));
        if numberbadpoints>loopimages
            numberbadpoints=loopimages
        end
        if numberbadpoints<1
            numberbadpoints=1
        end

%         displx(:,1)=-validx(:,1)+validx(:,numberbadpoints);
%         disply(:,1)=-validy(:,1)+validy(:,numberbadpoints);

        plot(validx(:,numberbadpoints),displx(:,numberbadpoints),'o');
        xlabel('position [pixel]')
        ylabel('displacement [pixel]')
        title(['Displacement versus position',sprintf(' (Current image #: %1g)',numberbadpoints)]);

%         validxtemp=validx;
%         validytemp=validy;
        displxtemp=displx;
        validxdelete=validxtemp;
        validydelete=validytemp;
        displxdelete=displxtemp;
        displydelete=displytemp;

        title(sprintf('Define the region of interest.  \n  All points ouside that region will be deleted'))

        [xgrid,ygrid]=ginput(2);
        x(1,1) = xgrid(1);
        x(1,2) = xgrid(2);
        y(1,1) = ygrid(2);
        y(1,2) = ygrid(1);

        deletepoints=find(validxdelete(:,numberbadpoints)>min(x) & validxdelete(:,numberbadpoints)<max(x) & displxdelete(:,numberbadpoints)<max(y) & displxdelete(:,numberbadpoints)>min(y));
        [loopnum one]=size(deletepoints);

        validxdelete(deletepoints,:)=[];
        validydelete(deletepoints,:)=[];


        % update temporary data matrices; delete bad points
        displxtemp(deletepoints,:)=[];
        displytemp(deletepoints,:)=[];
        validxtemp(deletepoints,:)=[];
        validytemp(deletepoints,:)=[];

        plot(validxtemp(:,numberbadpoints),displxtemp(:,numberbadpoints),'o');

        % delete point permanently?
        selectremove3 = menu(sprintf('Do you want to delete these markers permanently?'),'Yes','No');
        if selectremove3==1
            displx=displxtemp;
            disply=displytemp;
            validx=validxtemp;
            validy=validytemp;
        end
        if selectremove3==2
            displxtemp=displx;
            displytemp=disply;
            validxtemp=validx;
            validytemp=validy;
        end
        selectremove2 = menu(sprintf('Do you want to mark another bad point?'),'Yes','No');
        if selectremove2==2
            clear displx;
            clear disply;
            validxfirst=zeros(size(validx));
            validxfirst=validx(:,1)*ones(1,sizevalidx(1,2));
            validyfirst=zeros(size(validy));
            validyfirst=validy(:,1)*ones(1,sizevalidy(1,2));
            displx=validx-validxfirst;
            disply=validy-validyfirst;
            return
        end
    end

end

%---------------------------------
% Delete single points
% written by Chris
function [validx,validy]=removepoints_func3(validx,validy);

%sort out badpoints?

selection1 = menu(sprintf('Do you want to mark bad points?'),'Yes','No');
if selection1==2
    close all;
    return
end

% if yes
if selection1==1
    selection2=selection1;
    %     figure
    counter=0
    sizevalidx=size(validx);
    looppoints=sizevalidx(1,1);
    loopimages=sizevalidx(1,2)-1;
    defaultimage=loopimages
    numberbadpoints=0

    while selection2==1
        counter=counter+1
        clear xplot
        clear sizevalidx
        clear selection1
        clear selection2
        clear badpoints

        sizevalidx=size(validx);
        looppoints=sizevalidx(1,1);
        loopimages=sizevalidx(1,2)-1;

        clear displx;
        validxfirst=zeros(size(validx));
        validxfirst=validx(:,1)*ones(1,sizevalidx(1,2));
        displx=validx-validxfirst;

        % update temporary matrices
        displxtemp=displx;
        validxtemp=validx;
        validytemp=validy;
        %         resnormxtemp=resnormx;

        % get the image number from which the bad points will be chosen
        prompt = {'From which image do you want to choose the bad points?'};
        dlg_title = 'Bad points removal';
        num_lines= 1;
        if numberbadpoints==0
            defaultimage=loopimages
        end
        if numberbadpoints~0
            defaultimage=numberbadpoints
        end
        def     = {num2str(defaultimage)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        numberbadpoints = str2num(cell2mat(answer(1,1)));
        if numberbadpoints>loopimages
            numberbadpoints=loopimages
        end
        if numberbadpoints<1
            numberbadpoints=1
        end

        gridsizex=10*round(min(min(validx))/10):10:10*round(max(max(validx))/10);
        gridsizey=10*round(min(min(validy))/10):10:10*round(max(max(validy))/10);
        [XI,YI]=meshgrid(gridsizex,gridsizey);
        ZI=griddata(validx(:,numberbadpoints),validy(:,numberbadpoints),displx(:,numberbadpoints),XI,YI,'cubic');
        epsxx = gradient(ZI,10,10);

        % find max and min point and point them out
        mindisplx=find(displx(:,numberbadpoints)==min(displx(:,numberbadpoints)));
        maxdisplx=find(displx(:,numberbadpoints)==max(displx(:,numberbadpoints)));


        pcolor(XI,YI,epsxx);
        axis('equal')
        caxis([min(min(epsxx)) max(max(epsxx))])
        colorbar
        shading('interp')
        hold on
        plot3(validx(:,numberbadpoints),validy(:,numberbadpoints),displx(:,numberbadpoints)-min(displx(:,numberbadpoints)),'o','MarkerEdgeColor','k','MarkerFaceColor','g'), hold on;
        plot3(validx(mindisplx,numberbadpoints),validy(mindisplx,numberbadpoints),displx(mindisplx,numberbadpoints)-min(displx(:,numberbadpoints)),'o','MarkerEdgeColor','y','MarkerFaceColor','b')
        plot3(validx(maxdisplx,numberbadpoints),validy(maxdisplx,numberbadpoints),displx(maxdisplx,numberbadpoints)-min(displx(:,numberbadpoints)),'o','MarkerEdgeColor','y','MarkerFaceColor','r')
        axis([min(min(XI))-10 max(max(XI))+10 min(min(YI))-10 max(max(YI))+10])
        drawnow;
        hold off

        % get the bad point position

        title(sprintf('Click on the bad point',counter))
        [badpoint]=ginput(1);
        badpointx = badpoint(1,1);
        badpointy = badpoint(1,2);

        % find the point matching the given position
        wherethehellisthispoint=abs(validx(:,numberbadpoints)-badpoint(1,1))+abs(validy(:,numberbadpoints)-badpoint(1,2));
        badpointnum=find(wherethehellisthispoint==min(wherethehellisthispoint));

        % update temporary data matrices; delete bad points

        displxtemp(badpointnum,:)=[];
        validxtemp(badpointnum,:)=[];
        validytemp(badpointnum,:)=[];
        mindisplx=find(displxtemp(:,numberbadpoints)==min(displxtemp(:,numberbadpoints)));
        maxdisplx=find(displxtemp(:,numberbadpoints)==max(displxtemp(:,numberbadpoints)));

        % update the figure
        ZI=griddata(validxtemp(:,numberbadpoints),validytemp(:,numberbadpoints),displxtemp(:,numberbadpoints),XI,YI,'cubic');
        epsxx = gradient(ZI,10,10);
        pcolor(XI,YI,epsxx);
        axis('equal')
        caxis([min(min(epsxx)) max(max(epsxx))])
        colorbar
        shading('interp')
        hold on
        plot3(validxtemp(:,numberbadpoints),validytemp(:,numberbadpoints),displxtemp(:,numberbadpoints),'o','MarkerEdgeColor','k','MarkerFaceColor','g')
        plot3(validxtemp(mindisplx,numberbadpoints),validytemp(mindisplx,numberbadpoints),displxtemp(mindisplx,numberbadpoints)-min(displxtemp(:,numberbadpoints)),'o','MarkerEdgeColor','y','MarkerFaceColor','b')
        plot3(validxtemp(maxdisplx,numberbadpoints),validytemp(maxdisplx,numberbadpoints),displxtemp(maxdisplx,numberbadpoints)-min(displxtemp(:,numberbadpoints)),'o','MarkerEdgeColor','y','MarkerFaceColor','r')
        axis([min(min(XI))-10 max(max(XI))+10 min(min(YI))-10 max(max(YI))+10])
        drawnow;hold off;

        % delete point permanently?
        selection3 = menu(sprintf('Do you want to delete this point permanently?'),'Yes','No');
        if selection3==1
            displx=displxtemp;
            validx=validxtemp;
            validy=validytemp;
            %             resnormx=resnormxtemp;
        end
        if selection3==2
            displxtemp=displx;
            validxtemp=validx;
            validytemp=validy;
            %             resnormxtemp=resnormx;
        end
        selection2 = menu(sprintf('Do you want to mark another bad point?'),'Yes','No');
        if selection2==2
            close all;
            return
        end

    end
end

%---------------------------------
% Strain between two markers
% written by Chris
function [validx, validy,displx,disply] = strain_1D_2Points_func(validx, validy,displx,disply) ; % 1D strain calculation

sizevalidx=size(validx);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2)-1;
defaultimage=loopimages;
numberbadpoints=0;
clear selection3; selection3=1;

while selection3==1

    clear xplot
    clear sizevalidx
    clear selection1
    clear selection2
    clear badpoints

    sizevalidx=size(validx);
    looppoints=sizevalidx(1,1);
    loopimages=sizevalidx(1,2)-1;

    % update temporary matrices
    displxtemp=displx;
    validxtemp=validx;
    validytemp=validy;
    %     resnormxtemp=resnormx;

    % get the image number from which the bad points will be chosen
    prompt = {'Which image do you want for point selection?'};
    dlg_title = '1D Strain Plotting';
    num_lines= 1;
    if numberbadpoints==0
        defaultimage=loopimages;
    end
    if numberbadpoints~0
        defaultimage=numberbadpoints
    end
    def     = {num2str(defaultimage)};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    numberbadpoints = str2num(cell2mat(answer(1,1)));
    if numberbadpoints>loopimages
        numberbadpoints=loopimages;
    end
    if numberbadpoints<1
        numberbadpoints=1;
    end

    gridsizex=10*round(min(min(validx))/10):10:10*round(max(max(validx))/10);
    gridsizey=10*round(min(min(validy))/10):10:10*round(max(max(validy))/10);
    [XI,YI]=meshgrid(gridsizex,gridsizey);
    ZI=griddata(validx(:,numberbadpoints),validy(:,numberbadpoints),displx(:,numberbadpoints),XI,YI,'cubic');

    pcolor(XI,YI,ZI);
    axis('equal')
    caxis([min(min(ZI)) max(max(ZI))])
    colorbar
    shading('interp')
    hold on
    plot3(validx(:,numberbadpoints),validy(:,numberbadpoints),abs(displx(:,numberbadpoints)),'o','MarkerEdgeColor','k','MarkerFaceColor','g');
    axis([min(min(XI))-10 max(max(XI))+10 min(min(YI))-10 max(max(YI))+10])
    drawnow;

    % get the bad point position

    title(sprintf('Click on the two points for strain measurement'))
    [badpoint]=ginput(2);
    badpointx = badpoint(1,1);
    badpointy = badpoint(1,2);
    badpointx2 = badpoint(2,1);
    badpointy2 = badpoint(2,2);

    % find the point matching the given position
    wherethehellisthispoint=abs(validx(:,numberbadpoints)-badpoint(1,1))+abs(validy(:,numberbadpoints)-badpoint(1,2));
    badpointnum=find(wherethehellisthispoint==min(wherethehellisthispoint));
    wherethehellisthispoint2=abs(validx(:,numberbadpoints)-badpoint(2,1))+abs(validy(:,numberbadpoints)-badpoint(2,2));
    badpointnum2=find(wherethehellisthispoint2==min(wherethehellisthispoint2));


    % update the figure
    ZI=griddata(validxtemp(:,numberbadpoints),validytemp(:,numberbadpoints),displxtemp(:,numberbadpoints),XI,YI,'cubic');
    caxis([min(min(ZI)) max(max(ZI))])
    plot3(validxtemp(badpointnum,numberbadpoints),validytemp(badpointnum,numberbadpoints),displxtemp(badpointnum,numberbadpoints),'+','MarkerEdgeColor','k','MarkerFaceColor','g')
    plot3(validxtemp(badpointnum2,numberbadpoints),validytemp(badpointnum2,numberbadpoints),displxtemp(badpointnum2,numberbadpoints),'+','MarkerEdgeColor','k','MarkerFaceColor','r')
    hold off;
    axis([min(min(XI))-10 max(max(XI))+10 min(min(YI))-10 max(max(YI))+10])
    drawnow;

    epsilon1D=(displxtemp(badpointnum,:)-displxtemp(badpointnum2,:))/(validxtemp(badpointnum,1)-validxtemp(badpointnum2,1));
    epsilonsize=size(epsilon1D);
    figure; plot(1:epsilonsize(1,2),epsilon1D,'.');
    title(['True strain versus image from two picked points']);
    xlabel('Image number [ ]');
    ylabel('True Strain [ ]');


    selection3 = menu(sprintf('Do you want to choose 2 other points?'),'Yes','No');

    if selection3==2
        selection40 = menu(sprintf('Do you want to save the data as a text file?'),'Yes','No');
        if selection40==2
            return
        end
        if selection40==1
            numimagtemp = [1:epsilonsize(1,2)]';
            alltemp = [numimagtemp epsilon1D'];
            [FileNameBase,PathNameBase] = uiputfile('','Save file with image# vs. 1Dstrain');
            cd(PathNameBase)
            save(FileNameBase,'alltemp','-ASCII');
            %             save image_1Dstrain.txt alltemp -ASCII
            return
        end
    end
    close(gcf)
end
%---------------------------------
% Measure elastic slope
% written by Chris
function [validx, validy,displx,disply] = strain_1D_average_func(validx, validy,displx,disply) ; % 1D strain calculation
videoselection = menu(sprintf('Do you want to create a video?'),'Yes','No');
if videoselection==1
    mkdir('videostrain')
    cd('videostrain');
    Vid='Vid';
end
selection50=1;
validx_fit=validx;
displx_fit=displx;
minminvalidx=min(min(validx));
maxmaxvalidx=max(max(validx));
minminvalidy=min(min(validy));
maxmaxvalidy=max(max(validy));
minmindisplx=min(min(displx));
maxmaxdisplx=max(max(displx));
h= figure
while selection50==1
    %     figure
    [pointnumber imagenumber]=size(displx);
    for i=1:imagenumber;
        plot(validx_fit(:,i),displx_fit(:,i),'o');
        xdata=validx_fit(:,i);
        ydata=displx_fit(:,i);
        if i==1
            x(1)=0
            x(2)=0
        end
        [x,resnormx,residual,exitflagx,output]  = lsqcurvefit(@linearfit, [x(1) x(2)], xdata, ydata);
        hold on;
        ydatafit=x(1)*xdata+x(2);
        plot(xdata,ydatafit,'r');
        
        hold off
        slope(i,:)=[i x(1)];
        axis([minminvalidx maxmaxvalidx minmindisplx maxmaxdisplx])
        xlabel('position [pixel]')
        ylabel('displacement [pixel]')
        title(['Displacement versus position',sprintf(' (Current image #: %1g)',i)]);
        drawnow
        if videoselection==1
            u=i+10000;
            ustr=num2str(u);
            videoname=[Vid ustr '.jpg']
            saveas(h,videoname,'jpg')
        end
    end
    g1 = figure, plot(slope(:,1),slope(:,2));
    hold on
    plot(slope(:,1),slope(:,2),'.');
        xlabel('Image [ ]')
        ylabel('True Strain [ ]')
        title(['True Strain vs. Image #']);

    selection40 = menu(sprintf('Do you want to save the data as file?'),'Yes','No');
    if selection40==2

    end
    if selection40==1
        alltemp = [slope(:,1) slope(:,2)];
        [FileNameBase,PathNameBase] = uiputfile('','Save file with image# vs. 1Dstrain');
        cd(PathNameBase)
        save(FileNameBase,'alltemp','-ASCII');
        %         save image_1Dstrain_avg.txt alltemp -ASCII
    end

    selection50 = menu(sprintf('Do you want to analyse a selected area again?'),'Yes','No');
    if selection50==2
        clear validx_fit
        clear displx_fit
        return
    end
    if selection50==1
        close(g1)
        plot(validx_fit(:,imagenumber),displx_fit(:,imagenumber),'o');
        title(['True strain versus image from all markers']);
        xlabel('Image number [ ]');
        ylabel('True Strain [ ]');
        prompt = {'Min. x-position:','Max. x-position:'};
        dlg_title = 'Regime to be analyzed in pixels';
        num_lines= 1;
        def     = {'800','1200'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        minx= str2num(cell2mat(answer(1,1)));
        maxx= str2num(cell2mat(answer(2,1)));
        counter=0
        clear validx_fit
        clear displx_fit
        selectedmarkers=find(validx(:,imagenumber)>minx  & validx(:,imagenumber)<maxx);
        validx_fit=validx(selectedmarkers,:);
        displx_fit=displx(selectedmarkers,:);
        continue
    end

end


%---------------------------------
% 3D mesh plotting
% written by Chris
function [validx, validy,displx,disply]=meshplot(validx,validy,displx,disply);
h=figure;
sizevalidx=size(validx);
sizevalidy=size(validy);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);

videoselection = menu(sprintf('Do you want to create a video?'),'Yes','No');
if videoselection==1
    mkdir('video')
    cd('video');
    Vid='Vid';
end
gridsizex=10*round(min(min(validx))/10):10:10*round(max(max(validx))/10);
gridsizey=10*round(min(min(validy))/10):10:10*round(max(max(validy))/10);
[XI,YI]=meshgrid(gridsizex,gridsizey);
minminvalidx=min(min(validx));
maxmaxvalidx=max(max(validx));
minminvalidy=min(min(validy));
maxmaxvalidy=max(max(validy));
minmindisplx=min(min(displx));
maxmaxdisplx=max(max(displx));
minmindisply=min(min(disply));
maxmaxdisply=max(max(disply));
for i=1:(loopimages-1)
    ZI=griddata(validx(:,i),validy(:,i),displx(:,i),XI,YI,'cubic');
    mesh(XI,YI,ZI); hold on
    bottomplot=ones(size(validx))*minmindisplx;
    backyplaneplot=ones(size(validx))*maxmaxvalidy;
    backxplaneplot=ones(size(validx))*minminvalidx;
    plot3(validx(:,i),validy(:,i),displx(:,i),'.b');
    %         plot3(validx(:,i),backyplaneplot(:,i),displx(:,i),'.');
    plot3(backxplaneplot(:,i),validy(:,i),displx(:,i),'.g');
    xlabel('x-position [pixel]')
    ylabel('y-position [pixel]')
    zlabel('displacement [pixel]')
    %         plot3(validx(:,i),validy(:,i),bottomplot(:,i),'.');
    hold off
    title(['Displacement versus x-y-position',sprintf(' (Current image #: %1g)',i)]);
    axis([minminvalidx maxmaxvalidx minminvalidy maxmaxvalidy minmindisplx maxmaxdisplx])
    drawnow
    if videoselection==1
        u=i+10000;
        ustr=num2str(u);
        videoname=[Vid ustr '.jpg']
        saveas(h,videoname,'jpg')
    end
end

if videoselection==1
    cd('..')
end

%-------------------------------
% polyfit function
% written by Dan slightly changed by Chris
function [validx, validy,displx,disply]=polyfit3D(validx, validy,displx,disply);
close all
plot3dsurface_func(validx,validy,displx);

%---------------------------------------
% Just plot it
% written by Dan slightly changed by Chris
function plot3dsurface_func(validx,validy,displx,gridstyle,cropxx,cropyy)

sizevalidx=size(validx);
looppoints=sizevalidx(1,1);
loopimages=sizevalidx(1,2);
gridsizex=10*round(min(min(validx))/10):10:10*round(max(max(validx))/10);
gridsizey=10*round(min(min(validy))/10):10:10*round(max(max(validy))/10);
[XI,YI]=meshgrid(gridsizex,gridsizey);
ZI=griddata(validx(:,1),validy(:,1),displx(:,1),XI,YI,'cubic');
ZIsize=size(ZI);
displcolor = [-7 1];
straincolor = [-0.005 0.03];

maxminusminvalidx=(max(max(validx))-min(min(validx)));
maxminusminvalidy=(max(max(validx))-min(min(validy)));

for i=1:(loopimages-1)

    ZI=griddata(validx(:,i),validy(:,i),displx(:,i),XI,YI,'cubic');
    ZIsize=size(ZI);
    epsxx = gradient(ZI,(maxminusminvalidx/ZIsize(1,1)),(maxminusminvalidy/ZIsize(1,2)));

    subplot(2,1,1)
    pcolor(XI,YI,ZI)
    axis('equal')
    shading('interp')
    caxis(displcolor)
    h1 = colorbar;
    set(h1, 'PlotBoxAspectRatio',[2.0 10 8.0])
    set(h1, 'FontSize', 12);
    title(['Raw Displacement in x-direction',sprintf(' (Current image #: %1g)',i)]);

    subplot(2,1,2)
    pcolor(XI,YI,epsxx)
    axis('equal')
    shading('interp')
    caxis(straincolor)
    h1 = colorbar;
    set(h1, 'PlotBoxAspectRatio',[2.0 10 8.0])
    set(h1, 'FontSize', 12);
    title('Raw Strain in x-direction');

    drawnow

end

%--------------------------------------
% Delete some markers
% written by Chris
function [validx,validy] = removepoints_func2(validx,validy) ; %delete points

if exist('validx')==0
    [validx,Pathvalidx] = uigetfile('*.mat; *.txt','Open validx.mat or validx.txt');
    cd(Pathvalidx);
    validx=importdata(validx,'\t');
    [validy,Pathvalidy] = uigetfile('*.mat;*.txt','Open validy.mat or validy.txt');
    cd(Pathvalidy);
    validy=importdata(validy,'\t');
end


selectremove1 = menu(sprintf('Do you want to delete makers?'),'Yes','No');
if selectremove1==2

    return
end

% if yes
if selectremove1==1
    selectionremove2=selectremove1;
    %     figure
    counter=0
    sizevalidx=size(validx);
    looppoints=sizevalidx(1,1);
    loopimages=sizevalidx(1,2);
    defaultimage=loopimages
    numberbadpoints=0

    while selectionremove2==1
        counter=counter+1
        clear xplot
        clear sizevalidx
        clear selectremove11
        clear selection2
        %         clear badpoints

        sizevalidx=size(validx);
        looppoints=sizevalidx(1,1);
        loopimages=sizevalidx(1,2);

        % update temporary matrices
        %         displxtemp=displx;
        validxtemp=validx;
        validytemp=validy;
        %         resnormxtemp=resnormx;

        % get the image number from which the bad points will be chosen
        prompt = {'From which image do you want to delete markers?'};
        dlg_title = 'Marker removal';
        num_lines= 1;
        if numberbadpoints==0
            defaultimage=loopimages
        end
        if numberbadpoints~0
            defaultimage=numberbadpoints
        end
        def     = {num2str(defaultimage)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        numberbadpoints = str2num(cell2mat(answer(1,1)));
        if numberbadpoints>loopimages
            numberbadpoints=loopimages
        end
        if numberbadpoints<1
            numberbadpoints=1
        end

        displx(:,1)=-validx(:,1)+validx(:,numberbadpoints);
        displx(:,1)=displx(:,1)-min(displx(:,1));

        gridsizex=10*round(min(min(validx))/10):10:10*round(max(max(validx))/10);
        gridsizey=10*round(min(min(validy))/10):10:10*round(max(max(validy))/10);
        [XI,YI]=meshgrid(gridsizex,gridsizey);
        ZI=griddata(validx(:,numberbadpoints),validy(:,numberbadpoints),displx(:,1),XI,YI,'cubic');
        epsxx = gradient(ZI,10,10);

        pcolor(XI,YI,epsxx);
        axis('equal')
        caxis([min(min(epsxx)) max(max(epsxx))])
        colorbar
        shading('interp')
        hold on
        plot3(validx(:,numberbadpoints),validy(:,numberbadpoints),displx(:,1),'.','MarkerEdgeColor','k','MarkerFaceColor','g'), hold off;
        axis([min(min(XI))-10 max(max(XI))+10 min(min(YI))-10 max(max(YI))+10])
        drawnow;

        validxtemp=validx;
        validytemp=validy;
        displxtemp=displx;
        validxdelete=validxtemp;
        validydelete=validytemp;
        displxdelete=displxtemp;

        title(sprintf('Define the region of interest.  \n  All points ouside that region will be deleted'))

        [xgrid,ygrid]=ginput(2);
        x(1,1) = xgrid(1);
        x(1,2) = xgrid(2);
        y(1,1) = ygrid(2);
        y(1,2) = ygrid(1);

        deletepoints=find(validxdelete(:,numberbadpoints)>min(x) & validxdelete(:,numberbadpoints)<max(x) & validydelete(:,numberbadpoints)<max(y) & validydelete(:,numberbadpoints)>min(y));
        [loopnum one]=size(deletepoints);

        validxdelete(deletepoints,:)=[];
        validydelete(deletepoints,:)=[];

        plot3(validxtemp(:,numberbadpoints),validytemp(:,numberbadpoints),displxtemp(:,1),'o','MarkerEdgeColor','k','MarkerFaceColor','g'), hold off;

        % update temporary data matrices; delete bad points
        displxtemp(deletepoints,:)=[];
        validxtemp(deletepoints,:)=[];
        validytemp(deletepoints,:)=[];

        % update the figure
        gridsizex=10*round(min(min(validxtemp))/10):10:10*round(max(max(validxtemp))/10);
        gridsizey=10*round(min(min(validytemp))/10):10:10*round(max(max(validytemp))/10);
        [XI,YI]=meshgrid(gridsizex,gridsizey);
        ZI=griddata(validxtemp(:,numberbadpoints),validytemp(:,numberbadpoints),displxtemp(:,1),XI,YI,'cubic');
        epsxx = gradient(ZI,10,10);
        pcolor(XI,YI,epsxx);
        axis('equal')
        caxis([min(min(epsxx)) max(max(epsxx))])
        colorbar
        shading('interp')
        hold on
        plot3(validxtemp(:,numberbadpoints),validytemp(:,numberbadpoints),displxtemp(:,1),'o','MarkerEdgeColor','k','MarkerFaceColor','g'), hold off;
        axis([min(min(XI))-10 max(max(XI))+10 min(min(YI))-10 max(max(YI))+10])
        drawnow;

        % delete point permanently?
        selectremove3 = menu(sprintf('Do you want to delete these markers permanently?'),'Yes','No');
        if selectremove3==1
            displx=displxtemp;
            validx=validxtemp;
            validy=validytemp;
        end
        if selectremove3==2
            displxtemp=displx;
            validxtemp=validx;
            validytemp=validy;
        end
        selectremove2 = menu(sprintf('Do you want to mark another bad point?'),'Yes','No');
        if selectremove2==2
            clear displx;
            validxfirst=zeros(size(validx));
            validxfirst=validx(:,1)*ones(1,sizevalidx(1,2));
            displx=validx-validxfirst;
            return
        end
    end
end