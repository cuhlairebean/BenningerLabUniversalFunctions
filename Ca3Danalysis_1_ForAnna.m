% %% %%%% 3D zstack coorelation analysis %%%% %%
clearvars
clear Run_Nosilentcell
close all
clc

%%add universal functions
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
dbstop at 445 in Run_Nosilentcell

%% Set Path, Folder Names, Variables%%
tryloadingvalues = 1;%try loading masks - will try to load the data - boolean
cachannel = 1; %The number of the calcium channel
howmanychannel = 1; %put the number of channels total
zstacks = 8;

pcpath = 'E:\BarakBlum\'; %''
macpath = '/Volumes/Seagate Backup Plus Drive/BarakBlum/';
folder ={'KO 015\KO 015\693_UCN_3_Cre_ROBO_KO-015\'};

starttime = [-1]; %set to -1 if you want the first image to be the start otherwise put the number of the frame
endtime = [-1];%set to -1 if you want the last image

ending = '.xml';%file ending

resultsfolder = 'results';
savedappend = '';%adds extra writing to end of results file

if exist(pcpath, 'dir') == 7
    Path = pcpath;
    slash = '\';
    pc=1;
else
    Path = macpath;
    folder = regexprep(folder, '\', '/');
    slash = '/';
    pc=0;
end
time = datestr(datetime('today'));

istart = 1; %folder
jstart = 1; %file
kstart = 1; %islets

%% Run Analysis %%
errors = 0;
tic

for i = istart:length(folder)
    disp(['i = ' num2str(i)]);
    istart = 1;
    F = dir([Path folder{i} slash '*' ending]); %gets all the files within that folder with that ending
    
    for j=jstart:size(F,1)
        disp(['j = ' num2str(j)]);
        jstart = 1;
        name = F(j).name;
        name = regexprep(name, '[.].*', '');
        disp(name);
        fullpathname = [Path folder{i} slash name];
        filename.Location = [fullpathname ending];
        
        oLocation = pwd;
        cd([Path folder{i}])
        if (exist('results', 'dir')~=7)
            mkdir 'results';
        end
        cd(oLocation);
        x=0;
        
        if tryloadingvalues
             R = dir([Path folder{i} slash resultsfolder slash name '*.mat']);
            for k = kstart:size(R,1)
                kstart=1;
                try
                    load([Path folder{i} slash resultsfolder slash R(k).name]);
                catch
                    disp('Unable to load .mat file');
                end
            end
        end
        
        filename.Location = [fullpathname ending];
        results = [Path folder{i} slash resultsfolder slash name];
        
        
        [filename] = Run_3D_corrall(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i),ending);
        %[filename] = GetCorrelationMap(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i));
        
        savedfile = [results savedappend time];
        save([savedfile savedappend '.mat'], 'filename');
        saveAllFigsToPPT(savedfile)
        
        disp(['Analysis successful for ' fullpathname]);
        close all
        clear filename
    end
end
% %xlswrite([Path 'Phase-' datestr(datetime('today'))],difphase);
% save([Path 'corrall_workspace' datestr(datetime('today')) '.mat'])
% toc
