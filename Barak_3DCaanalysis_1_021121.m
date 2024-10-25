% %% %%%% 3D zstack coorelation analysis %%%% %%

clearvars
clear all
close all
clc


%%add universal functions
%dbstop at 185 in Run_3D_corrall

%% Set Path and Folder Names %%
tryloadingvalues = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%% set 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcpath ='G:\raw_data_intravital_2021\1_sec\';
macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/1_sec/';

% 1sec
% folder ={'32', '33'};%...
% cachannel = 2; %which channel is calcium
% howmanychannel = 3;
% zstacks = 1;
% starttime = [-1,-1];
% endtime = [-1,-1];
% ending = '.xml';

%% 1sec drift corrected
% folder ={'35_driftcorrect'};%...
% cachannel = 1; %which channel is calcium
% howmanychannel = 1;
% zstacks = 1;
% starttime = [-1];
% endtime = [-1];
% ending = '.tif';

%more 1sec
% folder ={'35', '36', '37', '38', '39', '40', '41', '42', '43'};%...
% cachannel = 3; %which channel is calcium
% howmanychannel = 3;
% zstacks = 1;
% starttime = [147,-1,-1,-1,143,190,69,-1,228];
% endtime = [375,175,-1,398,-1,-1,-1,-1,-1];
% ending = '.xml';

%more 6sec
% pcpath ='G:\raw_data_intravital_2021\6_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/6_sec/';
% 
% folder ={'20', '21', '22', '23', '24', '25', '26', '27', '29', '44'};%...
% cachannel = 3; %which channel is calcium
% howmanychannel = 3;
% zstacks = 3;
% starttime = [-1,-1,-1,-1,-1,-1,20,-1,42,-1];
% endtime =  [-1,-1,-1,-1,55,-1,-1,-1,-1,75];
% ending = '.xml';

%10sec
pcpath ='G:\raw_data_intravital_2021\10_sec\';
macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/10_sec/';
% 
% folder ={'1', '3', '4', '5', '6', '7', '8', '10', '11', '12', '13'};%...
% cachannel = 3; %which channel is calcium
% howmanychannel = 3;
% zstacks = 8;
% starttime = [-1,17,-1,-1,-1,-1,9,-1,-1,-1,-1];
% endtime =  [37,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
% ending = '.xml';
% 
folder ={'14', '15', '16','18', '19','611_27'};%...
cachannel = 3; %which channel is calcium
howmanychannel = 3;
zstacks = 8;
starttime = [30,-1,17,31,-1,21];
endtime =  [-1,-1,-1,-1,-1,-1];
ending = '.xml';

%%more 30sec
% pcpath ='G:\raw_data_intravital_2021\30_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/30_sec/';
% 
% folder ={ '45', '46', '47', '48'};%...
% cachannel = 3; %which channel is calcium
% howmanychannel = 3;
% zstacks = 12;
% starttime = [-1,-1,-1,-1];
% endtime =   [-1,18,-1,-1];
% ending = '.xml';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultsfolder = 'results';

if exist(pcpath, 'dir') == 7
    Path = pcpath;
    slash = '\';
    pc=1;
    addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
    folder = regexprep(folder, '/', '\');
else
    Path = macpath;
    folder = regexprep(folder, '\', '/');
    slash = '/';
    pc=0;
    addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/');
end
time = datestr(datetime('today'));

istart = 6; %folder
jstart = 1; %file

%% Run Analysis %%
errors = 0;
tic

%Txlsx = dir([Path folder '*Time.xlsx']);

%%create results folder if it does not exist

firstfile = 0;
funFreq = 0;

for i = istart:length(folder)
    disp(['i = ' num2str(i)]);
    istart = 1;
    F = dir([Path folder{i} slash '*' ending]);
    
    for j=jstart:size(F,1)
        disp(['j = ' num2str(j)]);
        jstart = 1;
        name = F(j).name;
        name = regexprep(name, '[.].*', '');
        disp(name);
        fullpathname = [Path folder{i} slash name];
        filename.Location = [fullpathname ending];
        %filename.TimeLoc = [Path folder Txlsx(1).name]
        results = [Path folder{i} slash resultsfolder slash name];
        filename.results = results;
        
        oLocation = pwd;
        cd([Path folder{i}])
        if (exist('results', 'dir')~=7)
            mkdir 'results';
        end
        cd(oLocation);
        x=0;
        
        R = dir([results '*2021.mat']);
        
        if tryloadingvalues
            for k = 1:size(R,1)
                try
                    analyzefile = contains(R(k).name,'stack','IgnoreCase',true);%doesnt analyze  files names snap
                    load([Path folder{i} slash resultsfolder slash R(k).name]);
                    filename.Location = [fullpathname ending];
                    filename.results = results;
                catch
                    disp('Unable to load .mat file');
                end
            end
        end
        filename.Location = [fullpathname ending];
        
         [filename] = Run_3D_corrall(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i),ending);
%        
        savedfile = [results '_corall' time];
        save([savedfile], 'filename', '-v7.3');
        saveAllFigsToPPT(savedfile)
        PhaseMap3D = filename.PhaseMap3D;
        FreqMap3D = filename.FreqMap3D;
        
        save([results 'PhaseMap.mat'], 'PhaseMap3D');
        save([results 'FreqMap.mat'], 'FreqMap3D');
%         
        [filename] = GetCorrelationMap(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i));
        
        disp(['Analysis successful for ' fullpathname]);
        close all
        clear filename
    end
end
% %xlswrite([Path 'Phase-' datestr(datetime('today'))],difphase);
% save([Path 'corrall_workspace' datestr(datetime('today')) '.mat'])
% toc
