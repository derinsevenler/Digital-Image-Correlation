function [displx disply]=large_displ;

% this function creats displx and disply which you need for
% automate_image_ld.m

%% choose resizing factor
% The reduction factor should be at least the largest step in your
% experiment divided by the corrsize you choose in cpcorr.m but will be
% better off being a little bit higher
reduction_factor=5;
prompt = {'Enter reduction factor - Image will be resized in the first run to track large displacement:'};
dlg_title = 'Reduction factor for large displacements';
num_lines= 1;
def     = {'5'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
reduction_factor = str2double(cell2mat(answer(1,1)));

%% choose a small grid for reduced size images
% the smaller the grid the faster this step
msgboxwicon=msgbox('Define grid for reduced image size - use 50 to 100 markers per image.')
waitfor(msgboxwicon)
[grid_x,grid_y]=grid_generator;
save grid_x_small.dat grid_x -ascii -tabs
save grid_y_small.dat grid_y -ascii -tabs

%% choose a larger grid for large sized images
msgboxwicon=msgbox('Define grid for detailed image analysis.')
waitfor(msgboxwicon)
[grid_x_full,grid_y_full]=grid_generator;

%% run automate_image.m
[validx,validy]=automate_image_rf(grid_x,grid_y,reduction_factor);

%% calculate displx and disply
disply=diff((mean(validy)-mean(validy(:,1)))*reduction_factor);
displx=diff((mean(validx)-mean(validx(:,1)))*reduction_factor);
displx=[0 displx];
disply=[0 disply];
save displx.dat displx -ascii -tabs
save disply.dat disply -ascii -tabs

%% start automate_image_ld
automate_image_ld(displx,disply,grid_x_full,grid_y_full);