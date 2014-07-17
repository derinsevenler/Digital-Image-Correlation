jobs=1;
jobnumber=0;

jobs = menu(sprintf('Do you want to add a job?'),...
    'Yes','No','Load','Cancel');

if jobs==3
    cd(matlabroot)
    load('job.mat');
    jobnumber=job{2,1};
    jobs = menu(sprintf('Do you want to add another job?'),...
        'Yes','No','Cancel');

    if jobs==3
        jobs==4;
    end
end



if jobs==1
    while jobs==1
        jobnumber=jobnumber+1;
        [FileNameBase,PathNameBase,filenamelist]=filelist_generator;
        job(1,jobnumber)={PathNameBase};
        job(2,1)={jobnumber};
        cd(PathNameBase);
        [grid_x,grid_y]=grid_generator(filenamelist(1,:),PathNameBase);
        jobs = menu(sprintf('Do you want to add another job?'),...
            'Yes','No','Save','Cancel');
        if jobs==3
            cd(matlabroot)
            save('job')
            jobs = menu(sprintf('Do you want to add another job?'),...
            'Yes','No','Cancel');
            if jobs==3
                jobs==4;
            end
        end

    end
end

if jobs==4
    return
end

if jobs==2;
    if jobnumber==0
        return
    end
    for i=1:jobnumber
        i
        job{1,i}
        cd(job{1,i})
        automate_image;
    end
end