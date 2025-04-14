% use to import all of the bdf files to letswave
clear all 

ERPpath = '/Users/florarosenzweig/Dropbox/Mon Mac (mac-114-318.local)/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/output/YOUNG/ERP';
logpath = fullfile(ERPpath, '/log');
outputpath = fullfile(ERPpath, '/Preprocessing');

%% Add letswave6 path
toolbox_path = '/Users/florarosenzweig/Dropbox/Mon Mac (mac-114-318.local)/Documents/MATLAB/Toolboxes/letswave6-master';
addpath(genpath(toolbox_path));

%% list TAP data folders and bdf files to import
cd (ERPpath)

% find all bdf files
list4import = dir('**/*sub0*ERP*.bdf');
condOrders = readcell(fullfile(outputpath,'condOrder.txt'));


%load('template_history.mat');

% load the bdf files in letswave and add condorder to name
for list = 1:size(list4import,1)
    cd(list4import(list).folder)
    disp(list4import(list).name)
    [out_header,out_data,message_string]= RLW_import_BDF(list4import(list).name);
    
    for subi = 2:size(condOrders,1) %start at 2 because 1 = Var1
        substring  = condOrders{subi,1};
        condstring = num2str(condOrders{subi,2});
        if contains(out_header.name,substring)
            out_headerCond = out_header;
            out_headerCond.name = strcat(out_header.name,'_',condstring);
            disp(['new file name: ', out_headerCond.name])
        end
    end
    
    %out_header.history = template_history;
    cd(outputpath)
    CLW_save([],out_headerCond,out_data);
    clear out_header out_headerCond out_data
end