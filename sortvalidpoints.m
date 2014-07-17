% written by Chris

function [validx,validy,resnormx]=sortvalidpoints(fitxy)

% This function will free the good peaks from fitxy, captured by
% peak_labelling, and send them back as validx and validy. The loops you
% will find are not absolutely necessary and reduce the speed of this
% function tremendously and will be changed probably in a future version if
% there is need.

sizefitxy=size(fitxy);
looppoints=sizefitxy(1,1);
loopimages=sizefitxy(1,2)/12;
validpoints=fitxy(:,sizefitxy(1,2)-11);
validnumber(:,:)=validpoints(1:(max(find(validpoints))));

for i=1:loopimages-1
    counter=1;
    for j=1:looppoints
        if counter<length(validnumber)+1
            if fitxy(j,i*12-11)==validnumber(counter)
                resnormx(counter,i)=fitxy(j,i*12);
                cropwidthx(counter,i)=fitxy(j,i*12-2);
                counter=counter+1;
            end
        end
    end
end

resnormxtemp(:,:)=resnormx(1:length(validnumber),:);
resnormx=resnormxtemp;
clear resnormxtemp;
normresxtest=mean((resnormx./cropwidthx)')';
number=1:1:looppoints;
figure, plot(1:length(validnumber),normresxtest,'.');
drawnow;

%set new normres
counter=1;
for i=1:length(validnumber)
    if abs(normresxtest(i,1))<10
        validnumbertest(counter,1)=validnumber(i,1);
        counter=counter+1;
    end
end

for i=1:loopimages-1
    counter=1;
    for j=1:looppoints
       if counter<length(validnumbertest)+1
            if fitxy(j,i*12-11)==validnumbertest(counter)
                validx(counter,i)=fitxy(j,i*12-9);
                validy(counter,i)=fitxy(j,i*12-5);
                resnormx(counter,i)=fitxy(j,i*12);
                counter=counter+1;
            end
        end
    end
end
pointnumber=size(validx);
validxtemp(:,:)=validx(1:pointnumber(1,1),:);
validx=validxtemp;
clear validxtemp;
validytemp(:,:)=validy(1:pointnumber(1,1),:);
validy=validytemp;
clear validytemp;
resnormxtemp(:,:)=resnormx(1:pointnumber(1,1),:);
clear resnormx;
resnormx=resnormxtemp;
